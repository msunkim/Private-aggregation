// Sigma-protocol spike.
// ---------------------------------------------------------------------------------------
// Generic non-interactive (Fiat--Shamir) proof of knowledge of a SHORT witness vector s
// satisfying a public linear relation  A s = u  over the ring  R_q = Z_q[x]/(x^N+1),
// with Lyubashevsky "rejection sampling" so the response leaks nothing about s.
//
// This is the common engine for all three PnF proofs:
//   POWF (sender well-formedness), PORK (aggregator re-randomizer), POCD (committee decryption);
// each is just a particular (A, u, s) instance of  A s = u  with bounded ||s||_inf.
//
// This spike validates the engine before instantiation:
//   [completeness]  an honest proof verifies;
//   [soundness]     a tampered response / wrong statement is rejected;
//   and prints proof size and prover/verifier time.
//
// NTL backend (Z_q[x]/(x^N+1)); the Fiat--Shamir hash here is a fast non-cryptographic
// stand-in for the random oracle (it does not affect the lattice-arithmetic timings we care
// about; a deployment would use a cryptographic hash).
//
// Build: see bench/CMakeLists.txt (ntl + gmp).

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include <iostream>
#include <vector>
#include <random>
#include <cstdint>
#include <chrono>
#include <cmath>

using namespace std;
using namespace NTL;
using Clock = std::chrono::steady_clock;
static double ms(Clock::time_point a, Clock::time_point b){ return chrono::duration<double,milli>(b-a).count(); }
static const char* PF(bool ok){ return ok ? "PASS" : "FAIL"; }

// ---------------------------------------------------------------- ring R_q = Z_q[x]/(x^N+1)
struct Ring {
    long N;
    ZZ q, half;
    ZZ_pXModulus F;                 // the modulus polynomial x^N + 1
    void init(long N_, const ZZ& q_) {
        N = N_; q = q_; half = q/2;
        ZZ_p::init(q);
        ZZ_pX f; SetCoeff(f, N, 1); SetCoeff(f, 0, 1);   // x^N + 1
        build(F, f);
    }
};

static inline long centered(const ZZ_p& c, const Ring& R) {
    ZZ v = rep(c);
    if (v > R.half) v -= R.q;
    return to_long(v);
}
static long infnorm(const ZZ_pX& a, const Ring& R) {
    long m = 0;
    for (long i = 0; i <= deg(a); ++i) { long v = labs(centered(coeff(a,i), R)); if (v > m) m = v; }
    return m;
}
static ZZ_pX unifPoly(long N) { ZZ_pX a; random(a, N); return a; }   // uniform in R_q
static ZZ_pX smallPoly(long N, long B, mt19937_64& rng) {            // coeffs uniform in [-B,B]
    ZZ_pX a; uniform_int_distribution<long> d(-B, B);
    for (long i = 0; i < N; ++i) SetCoeff(a, i, conv<ZZ_p>(conv<ZZ>(d(rng))));
    return a;
}

// ---------------------------------------------------------------- Fiat--Shamir challenge
// Hash the statement (u) and commitment (w) to a seed, then sample a challenge polynomial
// with exactly kappa nonzero coefficients, each +-1 (a sparse ternary challenge).
static uint64_t hashVec(const vector<ZZ_pX>& us, const vector<ZZ_pX>& ws, const Ring& R) {
    uint64_t h = 1469598103934665603ULL;             // FNV-1a (RO stand-in)
    auto absorb = [&](const vector<ZZ_pX>& vs){
        for (auto& v : vs) for (long i = 0; i < R.N; ++i) {
            uint64_t x = (uint64_t)centered(coeff(v,i), R);
            for (int b = 0; b < 8; ++b) { h ^= (x & 0xff); h *= 1099511628211ULL; x >>= 8; }
        }
    };
    absorb(us); absorb(ws);
    return h;
}
static ZZ_pX challenge(uint64_t seed, long N, long kappa) {
    mt19937_64 r(seed);
    ZZ_pX c;
    long placed = 0;
    while (placed < kappa) {
        long pos = (long)(r() % (uint64_t)N);
        if (!IsZero(coeff(c, pos))) continue;        // distinct positions
        SetCoeff(c, pos, conv<ZZ_p>(conv<ZZ>((r() & 1) ? 1 : -1)));
        ++placed;
    }
    return c;
}

// ---------------------------------------------------------------- linear-relation Sigma protocol
static vector<ZZ_pX> matVec(const vector<vector<ZZ_pX>>& A, const vector<ZZ_pX>& s, const Ring& R) {
    long m = A.size(); vector<ZZ_pX> out(m);
    for (long i = 0; i < m; ++i) { ZZ_pX acc; for (size_t j = 0; j < s.size(); ++j) { ZZ_pX t; MulMod(t, A[i][j], s[j], R.F); acc += t; } out[i] = acc; }
    return out;
}

struct Proof { ZZ_pX c; vector<ZZ_pX> z; long tries = 0; };

// Prove knowledge of short s (||s||_inf <= beta) with A s = u.  Bmask sets the mask range.
static Proof prove(const vector<vector<ZZ_pX>>& A, const vector<ZZ_pX>& u, const vector<ZZ_pX>& s,
                   const Ring& R, long beta, long kappa, long Bmask, mt19937_64& rng) {
    long d = s.size(), bound = Bmask - kappa*beta;
    Proof pf;
    while (true) {
        ++pf.tries;
        vector<ZZ_pX> y(d);
        for (long j = 0; j < d; ++j) y[j] = smallPoly(R.N, Bmask, rng);   // mask
        vector<ZZ_pX> w = matVec(A, y, R);
        ZZ_pX c = challenge(hashVec(u, w, R), R.N, kappa);
        vector<ZZ_pX> z(d); long zmax = 0;
        for (long j = 0; j < d; ++j) { ZZ_pX t; MulMod(t, c, s[j], R.F); z[j] = y[j] + t; long nn = infnorm(z[j], R); if (nn > zmax) zmax = nn; }
        if (zmax > bound) continue;                                       // reject -> restart
        pf.c = c; pf.z = z; return pf;
    }
}

static bool verify(const vector<vector<ZZ_pX>>& A, const vector<ZZ_pX>& u, const Proof& pf,
                   const Ring& R, long beta, long kappa, long Bmask) {
    long bound = Bmask - kappa*beta;
    for (auto& zj : pf.z) if (infnorm(zj, R) > bound) return false;        // norm check
    vector<ZZ_pX> Az = matVec(A, pf.z, R);
    long m = A.size(); vector<ZZ_pX> w(m);
    for (long i = 0; i < m; ++i) { ZZ_pX cu; MulMod(cu, pf.c, u[i], R.F); w[i] = Az[i] - cu; }
    return challenge(hashVec(u, w, R), R.N, kappa) == pf.c;                // recompute challenge
}

// proof size in bytes: z is d*N coefficients bounded by 'bound'; c is kappa sparse +-1 terms.
static double proofKB(const Proof& pf, const Ring& R, long beta, long kappa, long Bmask) {
    long bound = Bmask - kappa*beta;
    long zbits = (long)ceil(log2((double)(2*bound + 1)));
    double zbytes = (double)pf.z.size() * R.N * zbits / 8.0;
    double cbytes = (double)kappa * ((ceil(log2((double)R.N)) + 1) / 8.0);
    return (zbytes + cbytes) / 1024.0;
}

int main(int argc, char** argv) {
    const long N = (argc > 1) ? atol(argv[1]) : 1024;
    ZZ q; GenPrime(q, 59);
    Ring R; R.init(N, q);

    const long d = 3, m = 2;            // representative shape (e.g., A = [[b,1,0],[a,0,1]])
    const long beta = 20, kappa = 8;    // witness coeff bound; challenge weight
    const long Bmask = 2 * N * d * kappa * beta;

    mt19937_64 rng(20260622);
    std::cout << "=== Sigma-protocol spike (R_q = Z_q[x]/(x^" << N << "+1), log2 q = " << NumBits(q)
              << ") ===\n";
    std::cout << "shape: A is " << m << "x" << d << ", ||s||_inf <= " << beta
              << ", challenge weight kappa = " << kappa << "\n\n";

    // random instance: A uniform, s short, u = A s
    vector<vector<ZZ_pX>> A(m, vector<ZZ_pX>(d));
    for (long i = 0; i < m; ++i) for (long j = 0; j < d; ++j) A[i][j] = unifPoly(N);
    vector<ZZ_pX> s(d); for (long j = 0; j < d; ++j) s[j] = smallPoly(N, beta, rng);
    vector<ZZ_pX> u = matVec(A, s, R);

    auto t0 = Clock::now();
    Proof pf = prove(A, u, s, R, beta, kappa, Bmask, rng);
    auto t1 = Clock::now();
    bool ok = verify(A, u, pf, R, beta, kappa, Bmask);
    auto t2 = Clock::now();

    // soundness sanity: tamper the response, and try a wrong statement
    Proof bad = pf; bad.z[0] += conv<ZZ_pX>(1);
    bool rejTamper = !verify(A, u, bad, R, beta, kappa, Bmask);
    vector<ZZ_pX> uwrong = u; uwrong[0] += conv<ZZ_pX>(1);
    bool rejWrong = !verify(A, uwrong, pf, R, beta, kappa, Bmask);

    std::cout << "[completeness] honest proof verifies:     " << PF(ok) << "\n";
    std::cout << "[soundness]    tampered response rejected: " << PF(rejTamper) << "\n";
    std::cout << "[soundness]    wrong statement rejected:   " << PF(rejWrong) << "\n";
    std::cout << "\nprover " << ms(t0,t1) << " ms (" << pf.tries << " tries), verifier "
              << ms(t1,t2) << " ms, proof " << proofKB(pf, R, beta, kappa, Bmask) << " KB\n";

    bool all = ok && rejTamper && rejWrong;
    std::cout << "\n=== " << PF(all) << " ===\n";
    return all ? 0 : 1;
}
