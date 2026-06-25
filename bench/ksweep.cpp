// Committee-size (k) sweep for threshold BFV, fixed n = 1000 (depth 10).
// ---------------------------------------------------------------------------------------
// Measures the committee-dependent costs as the decryption committee size k grows:
//   setup     : k-party joint key generation + joint eval-mult-key (relinearization) chaining
//   p.d.      : partial decryptions (per member, and total over k members)
//   combine   : combining the k partial decryptions
// The homomorphic product, sender encryption, and recovery are k-independent and are not
// repeated here. Crucially, a BFV ciphertext is a single element at the modulus fixed by the
// depth (=10 for n=1000), regardless of how many factors were multiplied; so the committee
// costs below are identical to those for the full n=1000 aggregate, and we use a small product
// only to confirm correct multiset recovery (the eval-key chaining is right iff this passes).
//
// Generalizes the k=2 multiparty pattern (threshold_spike) to arbitrary k.
//
// Build: see bench/CMakeLists.txt (OpenFHE + ntl + gmp).

#include "openfhe.h"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

#include <iostream>
#include <vector>
#include <map>
#include <chrono>
#include <cmath>
#include <cstdint>

using namespace lbcrypto;
using Clock = std::chrono::steady_clock;
static double ms(Clock::time_point a, Clock::time_point b){ return std::chrono::duration<double,std::milli>(b-a).count(); }
static const char* PF(bool ok){ return ok ? "PASS" : "FAIL"; }

struct Committee {
    CryptoContext<DCRTPoly> cc;
    std::vector<KeyPair<DCRTPoly>> kp;     // k key pairs (chained)
    PublicKey<DCRTPoly> jointPk;
    double setupMs = 0;                     // k-party keygen + eval-mult-key chaining
};

// Build a k-party threshold context (depth fixed) and time the committee setup.
static Committee makeCommittee(int64_t t, uint32_t depth, uint32_t batch, int k) {
    CCParams<CryptoContextBFVRNS> p;
    p.SetPlaintextModulus(t);
    p.SetMultiplicativeDepth(depth);
    p.SetSecurityLevel(HEStd_128_classic);
    p.SetBatchSize(batch);
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(MULTIPARTY);

    Committee C; C.cc = cc;
    auto t0 = Clock::now();

    // ---- joint public key: chain MultipartyKeyGen over k parties ----
    C.kp.resize(k);
    C.kp[0] = cc->KeyGen();
    for (int i = 1; i < k; ++i) C.kp[i] = cc->MultipartyKeyGen(C.kp[i-1].publicKey);
    C.jointPk = C.kp[k-1].publicKey;
    auto tag = C.jointPk->GetKeyTag();

    // ---- joint eval-mult (relinearization) key ----
    // Step 1: combined key-switch key under the joint secret s = sum_i s_i.
    auto base = cc->KeySwitchGen(C.kp[0].secretKey, C.kp[0].secretKey);
    auto evalMultAB = base;
    for (int i = 1; i < k; ++i) {
        auto ek = cc->MultiKeySwitchGen(C.kp[i].secretKey, C.kp[i].secretKey, base);
        evalMultAB = cc->MultiAddEvalKeys(evalMultAB, ek, tag);
    }
    // Step 2: each party's MultiMultEvalKey, summed.
    auto evalMultFinal = cc->MultiMultEvalKey(C.kp[0].secretKey, evalMultAB, tag);
    for (int i = 1; i < k; ++i) {
        auto em = cc->MultiMultEvalKey(C.kp[i].secretKey, evalMultAB, tag);
        evalMultFinal = cc->MultiAddEvalMultKeys(evalMultFinal, em, tag);
    }
    cc->InsertEvalMultKey({evalMultFinal});

    C.setupMs = ms(t0, Clock::now());
    return C;
}

struct RunOut { double pdecTotal, pdecPer, fusion; bool ok; };

static RunOut runRound(const Committee& C, int64_t t, int k) {
    auto& cc = C.cc;
    const int nf = 16;                       // small product (correctness only; cost is n-independent)
    const long D = nf + 1;
    auto norm = [&](int64_t v){ return ((v % t) + t) % t; };

    std::vector<int64_t> roots(nf);
    for (int i = 0; i < nf; ++i) roots[i] = 1000 + i;   // distinct
    std::vector<Ciphertext<DCRTPoly>> cts;
    for (int i = 0; i < nf; ++i) {
        std::vector<int64_t> lanes(D);
        for (long j = 0; j < D; ++j) lanes[j] = norm((j + 1) - roots[i]);
        cts.push_back(cc->Encrypt(C.jointPk, cc->MakePackedPlaintext(lanes)));
    }
    while (cts.size() > 1) {
        std::vector<Ciphertext<DCRTPoly>> nxt;
        for (size_t i = 0; i + 1 < cts.size(); i += 2) nxt.push_back(cc->EvalMult(cts[i], cts[i+1]));
        if (cts.size() % 2 == 1) nxt.push_back(cts.back());
        cts = nxt;
    }
    auto prod = cts[0];

    // ---- partial decryptions: lead + (k-1) main ----
    std::vector<Ciphertext<DCRTPoly>> partials;
    auto a = Clock::now();
    auto lead = cc->MultipartyDecryptLead({prod}, C.kp[0].secretKey);
    partials.push_back(lead[0]);
    for (int i = 1; i < k; ++i) { auto m = cc->MultipartyDecryptMain({prod}, C.kp[i].secretKey); partials.push_back(m[0]); }
    auto b = Clock::now();

    // ---- fusion ----
    Plaintext res;
    cc->MultipartyDecryptFusion(partials, &res);
    res->SetLength(D);
    auto evals = res->GetPackedValue();
    for (auto& v : evals) v = norm(v);
    auto c = Clock::now();

    // ---- validate: interpolate + factor -> multiset matches ----
    NTL::ZZ_pPush push; NTL::ZZ_p::init(NTL::conv<NTL::ZZ>((long)t));
    NTL::vec_ZZ_p xa, ya; xa.SetLength(D); ya.SetLength(D);
    for (long j = 0; j < D; ++j) { xa[j]=NTL::conv<NTL::ZZ_p>((long)(j+1)); ya[j]=NTL::conv<NTL::ZZ_p>((long)evals[j]); }
    NTL::ZZ_pX M; NTL::interpolate(M, xa, ya); NTL::MakeMonic(M);
    NTL::vec_ZZ_p rts; NTL::FindRoots(rts, M);
    std::map<long,int> got, exp;
    for (long j = 0; j < rts.length(); ++j) got[NTL::conv<long>(NTL::rep(rts[j]))]++;
    for (int i = 0; i < nf; ++i) exp[(long)roots[i]]++;

    RunOut o;
    o.pdecTotal = ms(a,b);
    o.pdecPer   = ms(a,b) / k;
    o.fusion    = ms(b,c);
    o.ok        = (got == exp) && (rts.length() == nf);
    return o;
}

struct Stat{ double mean,sd; };
static Stat meanSD(const std::vector<double>& v){ double m=0; for(double x:v)m+=x; m/=v.size(); double s=0; for(double x:v)s+=(x-m)*(x-m); s=v.size()>1?std::sqrt(s/(v.size()-1)):0; return {m,s}; }

int main(int argc, char** argv) {
    const int TRIALS = (argc > 1) ? std::atoi(argv[1]) : 5;
    NTL::ZZ Mmod = NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start = ((NTL::conv<NTL::ZZ>((long)1) << 48) / Mmod) * Mmod + 1;
    while (!NTL::ProbPrime(start)) start += Mmod;
    const int64_t t = NTL::conv<long>(start);
    const uint32_t depth = 11, batch = 1024;          // n=1000 with +1 headroom (matches profile_b @ t~2^48)

    std::cout << "=== Committee-size sweep (threshold BFV, n=1000 => depth 11, N=2^15, "
              << TRIALS << " trials) ===\n\n";
    std::cout << "   k    setup (s)        p.d./member (ms)   p.d. total (ms)    combine (ms)   result\n";

    FILE* csv = std::fopen("results_ksweep.csv", "w");
    if (csv) std::fprintf(csv, "k,setup_mean_s,setup_sd_s,pdec_per_mean_ms,pdec_per_sd_ms,pdec_total_mean_ms,pdec_total_sd_ms,combine_mean_ms,combine_sd_ms,ok\n");

    for (int k : {2, 5, 10, 20}) {
        std::vector<double> su, pp, pt, fu; bool ok = true;
        for (int tr = 0; tr < TRIALS; ++tr) {
            Committee C = makeCommittee(t, depth, batch, k);
            RunOut o = runRound(C, t, k);
            su.push_back(C.setupMs/1000.0); pp.push_back(o.pdecPer); pt.push_back(o.pdecTotal); fu.push_back(o.fusion); ok &= o.ok;
        }
        Stat ss=meanSD(su), sp=meanSD(pp), st=meanSD(pt), sf=meanSD(fu);
        printf("%4d   %6.2f+/-%-5.2f   %7.2f+/-%-6.2f   %8.1f+/-%-6.1f   %6.2f+/-%-5.2f   %s\n",
               k, ss.mean,ss.sd, sp.mean,sp.sd, st.mean,st.sd, sf.mean,sf.sd, PF(ok));
        fflush(stdout);
        if (csv) { std::fprintf(csv, "%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n",
                    k, ss.mean,ss.sd, sp.mean,sp.sd, st.mean,st.sd, sf.mean,sf.sd, ok?1:0); std::fflush(csv); }
    }
    if (csv) std::fclose(csv);
    std::cout << "\n=== done (archive: results_ksweep.csv) ===\n";
    return 0;
}
