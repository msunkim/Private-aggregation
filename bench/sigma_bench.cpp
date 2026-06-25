// Sigma-protocol benchmark: POWF / PORK / POCD.
// ---------------------------------------------------------------------------------------
// Builds on the engine validated in sigma_spike: a Fiat--Shamir proof of knowledge of a
// short witness s with A s = u over R_q = Z_q[x]/(x^N+1), Lyubashevsky rejection sampling.
// Here we (i) use PER-COORDINATE norm bounds, since each relation mixes small witnesses
// (encryption randomness u, errors e0,e1, key share s_j, all ||.||<=B_err) with a few
// LARGE ones, and (ii) use a large proof modulus q_zk (bignum) to avoid wraparound.
//
// Representative instantiations (one large prime q_zk; matching the FHE ring dimension N):
//   POWF (sender):     c0 = b u + e0 + Delta*m ,  c1 = a u + e1
//                      witness (u, e0, e1, m); m is the encoded message, bounded by t ~ 2^48
//                      (a bounded-message coordinate), the rest small.
//   PORK (aggregator): difference of products is a *sanitizing* re-randomization of 0
//                      witness (u', e0', e1'); sanitization floods the error, so e0',e1'
//                      are bounded by B_fl ~ 2^40 (no longer short), u' small.
//   POCD (committee):  d_j = c1 * s_j + e_smudge
//                      witness (s_j small, e_smudge), where e_smudge is the threshold
//                      flooding noise, bounded by the statistical smudging gap ~ 2^40.
//
// Reports prover time, verifier time, #rejection tries, and proof size for each, as
// mean +- sample SD over several trials. The aggregator additionally performs n POWF
// verifications per round; we report the single-proof cost and note the n-scaling.
//
// Build: see bench/CMakeLists.txt (ntl + gmp).

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>
#include <cstdint>

using namespace std;
using namespace NTL;
using Clock = std::chrono::steady_clock;
static double ms(Clock::time_point a, Clock::time_point b){ return chrono::duration<double,milli>(b-a).count(); }

// ---------------------------------------------------------------- ring R_q = Z_q[x]/(x^N+1)
struct Ring {
    long N; ZZ q, half; ZZ_pXModulus F;
    void init(long N_, const ZZ& q_) {
        N = N_; q = q_; half = q/2;
        ZZ_p::init(q);
        ZZ_pX f; SetCoeff(f, N, 1); SetCoeff(f, 0, 1);
        build(F, f);
    }
};
static inline ZZ centered(const ZZ_p& c, const Ring& R){ ZZ v = rep(c); if (v > R.half) v -= R.q; return v; }
static ZZ infnorm(const ZZ_pX& a, const Ring& R){ ZZ m; m=0; for (long i=0;i<=deg(a);++i){ ZZ v=abs(centered(coeff(a,i),R)); if (v>m) m=v; } return m; }
static ZZ_pX unifPoly(long N){ ZZ_pX a; random(a,N); return a; }
static ZZ_pX smallPoly(long N, const ZZ& B){                       // coeffs uniform in [-B,B]
    ZZ_pX a; ZZ span = 2*B + 1;
    for (long i=0;i<N;++i){ ZZ v = RandomBnd(span) - B; SetCoeff(a,i, conv<ZZ_p>(v)); }
    return a;
}

// ---------------------------------------------------------------- Fiat--Shamir challenge
static uint64_t hashVec(const vector<ZZ_pX>& us, const vector<ZZ_pX>& ws, const Ring& R){
    uint64_t h = 1469598103934665603ULL; unsigned char buf[16];
    auto absorb = [&](const vector<ZZ_pX>& vs){
        for (auto& v : vs) for (long i=0;i<R.N;++i){
            BytesFromZZ(buf, rep(coeff(v,i)), 12);
            for (int b=0;b<12;++b){ h ^= buf[b]; h *= 1099511628211ULL; }
        }
    };
    absorb(us); absorb(ws); return h;
}
static ZZ_pX challenge(uint64_t seed, long N, long kappa){
    // deterministic RO stand-in: splitmix64 stream -> kappa distinct +-1 positions
    ZZ_pX c; long placed=0; uint64_t s=seed;
    auto next=[&](){ s+=0x9E3779B97F4A7C15ULL; uint64_t z=s; z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31); };
    while (placed<kappa){ long pos=(long)(next()%(uint64_t)N); if (!IsZero(coeff(c,pos))) continue; SetCoeff(c,pos, conv<ZZ_p>(ZZ((next()&1)?1:-1))); ++placed; }
    return c;
}

// ---------------------------------------------------------------- linear-relation Sigma protocol
static vector<ZZ_pX> matVec(const vector<vector<ZZ_pX>>& A, const vector<ZZ_pX>& s, const Ring& R){
    long m=A.size(); vector<ZZ_pX> out(m);
    for (long i=0;i<m;++i){ ZZ_pX acc; for (size_t j=0;j<s.size();++j){ ZZ_pX t; MulMod(t,A[i][j],s[j],R.F); acc+=t; } out[i]=acc; }
    return out;
}
struct Proof { ZZ_pX c; vector<ZZ_pX> z; vector<ZZ_pX> w; long tries=0; };  // w = commitment (for batch verify)

struct Relation {
    string name;
    vector<vector<ZZ_pX>> A;   // m x d
    vector<ZZ_pX> u, s;        // m ; d
    vector<ZZ> beta, Bmask;    // per-coordinate witness bound / mask range (d)
    long kappa;
};

static Proof prove(const Relation& rel, const Ring& R){
    long d=rel.s.size(); Proof pf;
    while (true){
        ++pf.tries;
        vector<ZZ_pX> y(d); for (long j=0;j<d;++j) y[j]=smallPoly(R.N, rel.Bmask[j]);
        vector<ZZ_pX> w = matVec(rel.A, y, R);
        ZZ_pX c = challenge(hashVec(rel.u, w, R), R.N, rel.kappa);
        vector<ZZ_pX> z(d); bool ok=true;
        for (long j=0;j<d;++j){ ZZ_pX t; MulMod(t,c,rel.s[j],R.F); z[j]=y[j]+t; if (infnorm(z[j],R) > rel.Bmask[j]-rel.kappa*rel.beta[j]){ ok=false; break; } }
        if (!ok) continue;
        pf.c=c; pf.z=z; pf.w=w; return pf;
    }
}
static bool verify(const Relation& rel, const Proof& pf, const Ring& R){
    long d=rel.s.size();
    for (long j=0;j<d;++j) if (infnorm(pf.z[j],R) > rel.Bmask[j]-rel.kappa*rel.beta[j]) return false;
    vector<ZZ_pX> Az = matVec(rel.A, pf.z, R);
    long m=rel.A.size(); vector<ZZ_pX> w(m);
    for (long i=0;i<m;++i){ ZZ_pX cu; MulMod(cu,pf.c,rel.u[i],R.F); w[i]=Az[i]-cu; }
    return challenge(hashVec(rel.u, w, R), R.N, rel.kappa) == pf.c;
}
static double proofKB(const Relation& rel, const Ring& R){
    double bits=0;
    for (long j=0;j<(long)rel.s.size();++j){ ZZ bound = rel.Bmask[j]-rel.kappa*rel.beta[j]; bits += (double)R.N * NumBits(2*bound+1); }
    bits += rel.kappa * (ceil(log2((double)R.N)) + 1);
    return bits/8.0/1024.0;
}

// ---------------------------------------------------------------- relation builders
static const long C_MASK = 10;     // mask factor: B = C_MASK * N * kappa * beta  (~0.9 accept/coord)
static Relation makeRel(const string& nm, long m, long d, const vector<ZZ>& beta, long kappa, const Ring& R){
    Relation rel; rel.name=nm; rel.kappa=kappa; rel.beta=beta; rel.Bmask.resize(d);
    for (long j=0;j<d;++j) rel.Bmask[j] = ZZ(C_MASK) * R.N * kappa * beta[j];
    rel.A.assign(m, vector<ZZ_pX>(d));
    for (long i=0;i<m;++i) for (long j=0;j<d;++j) rel.A[i][j]=unifPoly(R.N);   // public (pk-derived) entries
    rel.s.resize(d); for (long j=0;j<d;++j) rel.s[j]=smallPoly(R.N, beta[j]);  // short witness
    rel.u = matVec(rel.A, rel.s, R);
    return rel;
}

// Batch-verify n proofs sharing the SAME matrix A but distinct statements u_i. Each proof
// carries its commitment w_i. Cheap per-proof checks (norms + c_i = H(u_i,w_i)) plus ONE
// batched algebraic check  A(sum_i rho_i z_i) == sum_i rho_i (w_i + c_i u_i)  with random
// scalars rho_i in Z_q (soundness error ~ 1/q). Replaces n matrix-vector products by one.
static bool batchVerify(const vector<vector<ZZ_pX>>& A, const vector<vector<ZZ_pX>>& us,
                        const vector<Proof>& pfs, const Ring& R,
                        const vector<ZZ>& beta, const vector<ZZ>& Bmask, long kappa) {
    long n=pfs.size(), d=A[0].size(), m=A.size();
    for (long i=0;i<n;++i){
        for (long j=0;j<d;++j) if (infnorm(pfs[i].z[j],R) > Bmask[j]-kappa*beta[j]) return false;
        if (challenge(hashVec(us[i], pfs[i].w, R), R.N, kappa) != pfs[i].c) return false;
    }
    vector<ZZ_pX> zc(d), rhs(m);
    for (long i=0;i<n;++i){
        ZZ_p rho; random(rho);
        for (long j=0;j<d;++j){ ZZ_pX t; mul(t, pfs[i].z[j], rho); zc[j]+=t; }
        for (long r=0;r<m;++r){ ZZ_pX cu; MulMod(cu, pfs[i].c, us[i][r], R.F); ZZ_pX t; mul(t, pfs[i].w[r]+cu, rho); rhs[r]+=t; }
    }
    vector<ZZ_pX> lhs = matVec(A, zc, R);
    for (long r=0;r<m;++r) if (lhs[r]!=rhs[r]) return false;
    return true;
}

struct Stat{ double mean,sd; };
static Stat meanSD(const vector<double>& v){ double m=0; for(double x:v)m+=x; m/=v.size(); double s=0; for(double x:v)s+=(x-m)*(x-m); s=v.size()>1?sqrt(s/(v.size()-1)):0; return {m,s}; }

int main(int argc, char** argv){
    const long N = (argc>1)?atol(argv[1]):32768;
    const int TRIALS = (argc>2)?atoi(argv[2]):5;
    ZZ q; GenPrime(q, 80);                       // proof modulus q_zk
    Ring R; R.init(N, q);
    SetSeed(conv<ZZ>(20260622));

    const long kappa = 8;
    const ZZ Berr  = ZZ(20);                      // small error/ternary bound
    const ZZ Bmsg  = power2_ZZ(48);               // message m ~ t ~ 2^48 (root = <=8-bit category || 40-bit tag)
    const ZZ Bsmg  = power2_ZZ(40);               // smudging (flooding) noise ~ 2^40 (lambda=40)
    const ZZ Bfl   = power2_ZZ(40);               // ciphertext sanitization flooding on c^x ~ 2^40 (lambda=40)

    cout << "=== Sigma-protocol benchmark (R_q = Z_q[x]/(x^" << N << "+1), log2 q_zk = " << NumBits(q)
         << ", kappa = " << kappa << ", " << TRIALS << " trials) ===\n\n";
    cout << "relation   shape   prover (ms)        verifier (ms)      tries   proof (KB)\n";

    struct Spec{ string nm; long m,d; vector<ZZ> beta; };
    vector<Spec> specs = {
        {"POWF", 2, 4, {Berr,Berr,Berr,Bmsg}},   // (u,e0,e1,m) ; m bounded by t
        {"PORK", 2, 3, {Berr,Bfl,Bfl}},          // (u',e0',e1'): sanitizing re-randomizer, errors flooded to B_fl
        {"POCD", 1, 2, {Berr,Bsmg}},             // (s_j, e_smudge)
    };

    for (auto& sp : specs){
        Relation rel = makeRel(sp.nm, sp.m, sp.d, sp.beta, kappa, R);
        vector<double> tp, tv; vector<double> tr; bool ok=true;
        double kb = proofKB(rel, R);
        for (int t=0;t<TRIALS;++t){
            auto a=Clock::now(); Proof pf=prove(rel,R); auto b=Clock::now();
            bool v=verify(rel,pf,R); auto c=Clock::now();
            tp.push_back(ms(a,b)); tv.push_back(ms(b,c)); tr.push_back(pf.tries); ok &= v;
        }
        Stat sp_=meanSD(tp), sv=meanSD(tv), st=meanSD(tr);
        printf("%-6s   %ldx%ld   %7.1f+/-%-6.1f  %7.1f+/-%-6.1f  %4.1f   %8.1f   %s\n",
               sp.nm.c_str(), sp.m, sp.d, sp_.mean,sp_.sd, sv.mean,sv.sd, st.mean, kb, ok?"ok":"FAIL");
    }
    // ---- batch verification: aggregator verifies n POWF proofs (same A, distinct u_i) ----
    {
        Relation rel = makeRel("POWF", 2, 4, {Berr,Berr,Berr,Bmsg}, kappa, R);
        Proof pf = prove(rel, R);
        auto a=Clock::now(); bool sv=verify(rel, pf, R); auto b=Clock::now();
        double tsingle = ms(a,b);
        cout << "\n=== Batch verification (POWF, shared A): n single verifications vs one batched ===\n";
        cout << "    n   n x single (ms)   batched (ms)   speedup   ok\n";
        for (long n : {50L, 300L, 1000L}) {
            vector<vector<ZZ_pX>> us(n, rel.u);     // identical valid instances (timing-representative)
            vector<Proof> pfs(n, pf);
            auto c=Clock::now(); bool bok = batchVerify(rel.A, us, pfs, R, rel.beta, rel.Bmask, kappa); auto d=Clock::now();
            double tb = ms(c,d), tsn = tsingle * n;
            printf("  %4ld   %12.1f   %11.1f   %6.1fx   %s\n", n, tsn, tb, tsn/tb, (sv&&bok)?"ok":"FAIL");
        }
        // soundness sanity: one tampered proof must make the batch reject
        vector<vector<ZZ_pX>> us(8, rel.u); vector<Proof> pfs(8, pf);
        pfs[3].z[0] += conv<ZZ_pX>(1);
        bool rej = !batchVerify(rel.A, us, pfs, R, rel.beta, rel.Bmask, kappa);
        cout << "  [soundness] batch rejects a tampered proof: " << (rej?"PASS":"FAIL") << "\n";
    }
    return 0;
}
