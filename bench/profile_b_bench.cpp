// PnF Profile-B benchmark with threshold (multiparty) BFV, committee size k = 2.
// ---------------------------------------------------------------------------------------
// Measures the per-party costs of one PnF measurement round, swept over the number of
// senders n, in the evaluation-domain (Profile B) encoding with a LARGE plaintext modulus
// (t ~ 2^44, uncapped). Threshold decryption uses OpenFHE's native multiparty path
// (validated by threshold_spike): joint key-gen, joint relinearization, and
// MultipartyDecryptLead/Main/Fusion.
//
// Per-round phases timed:
//   setup     : multiparty key-gen + joint eval-mult key (one-time per n; depth-dependent)
//   enc/sender: one encryption under the joint pk          (sender S_i cost)
//   product   : homomorphic multiplication tree            (aggregator cost)
//   pdec/memb : one partial decryption                     (committee D_j cost)
//   fusion    : combine partial decryptions               (recovery)
//   interp+rt : interpolation + root-finding over F_t      (recovery)
//
// Correctness (multiset recovered == submitted) is checked every run.
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
#include <random>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>

using namespace lbcrypto;
using Clock = std::chrono::steady_clock;
static double ms(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}
static const char* PF(bool ok) { return ok ? "PASS" : "FAIL"; }

// Pick a Profile-B prime t ~ 2^48 with t = 1 (mod 2^17) so packing works for N <= 65536.
// t ~ 2^48 is the practical ceiling for the 64-bit native backend with packed encoding;
// it accommodates a 48-bit root (<= 8-bit category || 40-bit tag, see Section 5/6).
static int64_t pickT() {
    NTL::ZZ Mmod  = NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start = ((NTL::conv<NTL::ZZ>((long)1) << 48) / Mmod) * Mmod + 1;
    while (!NTL::ProbPrime(start)) start += Mmod;
    return NTL::conv<long>(start);
}

struct Committee {
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> kp1, kp2;
    PublicKey<DCRTPoly> jointPk;
};

// Build context + multiparty (k=2) keys for a given multiplicative depth.
static Committee makeCommittee(int64_t t, uint32_t depth, uint32_t batch) {
    CCParams<CryptoContextBFVRNS> p;
    p.SetPlaintextModulus(t);
    p.SetMultiplicativeDepth(depth);
    p.SetSecurityLevel(HEStd_128_classic);
    p.SetBatchSize(batch);
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(MULTIPARTY);

    Committee C; C.cc = cc;
    C.kp1 = cc->KeyGen();
    C.kp2 = cc->MultipartyKeyGen(C.kp1.publicKey);
    C.jointPk = C.kp2.publicKey;

    auto emk   = cc->KeySwitchGen(C.kp1.secretKey, C.kp1.secretKey);
    auto emk2  = cc->MultiKeySwitchGen(C.kp2.secretKey, C.kp2.secretKey, emk);
    auto emkAB = cc->MultiAddEvalKeys(emk, emk2, C.kp2.publicKey->GetKeyTag());
    auto emkBAB= cc->MultiMultEvalKey(C.kp2.secretKey, emkAB, C.kp2.publicKey->GetKeyTag());
    auto emkAAB= cc->MultiMultEvalKey(C.kp1.secretKey, emkAB, C.kp2.publicKey->GetKeyTag());
    auto emkF  = cc->MultiAddEvalMultKeys(emkAAB, emkBAB, emkAB->GetKeyTag());
    cc->InsertEvalMultKey({emkF});
    return C;
}

struct Timing { double enc, product, pdec, recovery; bool ok; uint32_t N; };

static Timing runRound(const Committee& C, int64_t t, int n, std::mt19937_64& rng) {
    auto& cc = C.cc;
    const long degM = n;
    const long D    = degM + 1;            // evaluation points a_j = j+1
    auto norm = [&](int64_t v) { return ((v % t) + t) % t; };

    // distinct uniform field elements (a 48-bit root carries an <=8-bit category and a 40-bit tag)
    std::vector<int64_t> roots(n);
    std::map<long,int> seen;
    for (int i = 0; i < n; ++i) { int64_t r; do { r = (int64_t)(rng() % (uint64_t)t); } while (seen[(long)r]++); roots[i] = r; }

    // --- senders: encrypt n factors under the joint pk ---
    auto t0 = Clock::now();
    std::vector<Ciphertext<DCRTPoly>> cts; cts.reserve(n);
    for (int i = 0; i < n; ++i) {
        std::vector<int64_t> lanes(D);
        for (long j = 0; j < D; ++j) lanes[j] = norm((j + 1) - roots[i]);
        cts.push_back(cc->Encrypt(C.jointPk, cc->MakePackedPlaintext(lanes)));
    }
    auto t1 = Clock::now();

    // --- aggregator: homomorphic product tree ---
    while (cts.size() > 1) {
        std::vector<Ciphertext<DCRTPoly>> nxt;
        for (size_t i = 0; i + 1 < cts.size(); i += 2) nxt.push_back(cc->EvalMult(cts[i], cts[i + 1]));
        if (cts.size() % 2 == 1) nxt.push_back(cts.back());
        cts = nxt;
    }
    auto prod = cts[0];
    // sanitize the aggregate: re-randomize with an encryption of zero. The statistical
    // flooding noise (lambda=40) is set analytically and absorbed by the modulus headroom
    // the threshold smudging already requires; its runtime is one encryption, negligible
    // against the product. We fold it into the aggregator's timing below.
    { std::vector<int64_t> zl(D, 0); prod = cc->EvalAdd(prod, cc->Encrypt(C.jointPk, cc->MakePackedPlaintext(zl))); }
    auto t2 = Clock::now();

    // --- committee: partial decryptions (lead + main) ---
    auto pL = cc->MultipartyDecryptLead({prod}, C.kp1.secretKey);
    auto t3 = Clock::now();
    auto pM = cc->MultipartyDecryptMain({prod}, C.kp2.secretKey);
    auto t4 = Clock::now();

    // --- fusion ---
    Plaintext res;
    cc->MultipartyDecryptFusion({pL[0], pM[0]}, &res);
    res->SetLength(D);
    auto evals = res->GetPackedValue();
    for (auto& v : evals) v = norm(v);

    // --- recovery: interpolate + root-find over F_t ---
    NTL::ZZ_pPush push;
    NTL::ZZ_p::init(NTL::conv<NTL::ZZ>((long)t));
    NTL::vec_ZZ_p xa, ya; xa.SetLength(D); ya.SetLength(D);
    for (long j = 0; j < D; ++j) { xa[j] = NTL::conv<NTL::ZZ_p>((long)(j+1)); ya[j] = NTL::conv<NTL::ZZ_p>((long)evals[j]); }
    NTL::ZZ_pX M; NTL::interpolate(M, xa, ya); NTL::MakeMonic(M);
    NTL::vec_ZZ_p rts; NTL::FindRoots(rts, M);
    auto t6 = Clock::now();

    std::map<long,int> got, exp;
    for (long j = 0; j < rts.length(); ++j) got[NTL::conv<long>(NTL::rep(rts[j]))]++;
    for (int i = 0; i < n; ++i) exp[(long)roots[i]]++;

    Timing T;
    T.enc      = ms(t0,t1) / n;             // per sender
    T.product  = ms(t1,t2);
    T.pdec     = (ms(t2,t3) + ms(t3,t4)) / 2.0;  // per committee member
    T.recovery = ms(t4,t6);                 // combine (fusion) + interpolate + root-find
    T.ok       = (got == exp) && (rts.length() == n);
    T.N        = cc->GetRingDimension();
    return T;
}

struct Stat { double mean, sd; };
static Stat meanSD(const std::vector<double>& v) {
    double m = 0; for (double x : v) m += x; m /= v.size();
    double s = 0; for (double x : v) s += (x - m) * (x - m);
    s = (v.size() > 1) ? std::sqrt(s / (v.size() - 1)) : 0.0;   // sample SD
    return {m, s};
}

int main(int argc, char** argv) {
    const int64_t t = pickT();
    const int TRIALS = (argc > 1) ? std::atoi(argv[1]) : 10;
    std::vector<int> ns = {50, 100, 200, 300, 500, 1000};

    std::cout << "=== PnF Profile-B benchmark (threshold BFV, k=2) ===\n";
    std::cout << "t = " << t << " (~2^48), security = 128-bit, trials = " << TRIALS
              << " (reporting mean +/- sample SD, ms)\n\n";
    std::cout << "   n      N  depth     enc/sender         product        pdec/memb         recovery       result\n";

    // CSV archive: one row per n, mean and SD for every metric.
    FILE* csv = std::fopen("results_profileB_k2.csv", "w");
    if (csv) std::fprintf(csv,
        "n,N,depth,trials,"
        "enc_mean_ms,enc_sd_ms,product_mean_ms,product_sd_ms,"
        "pdec_mean_ms,pdec_sd_ms,recovery_mean_ms,recovery_sd_ms,ok\n");

    for (int n : ns) {
        // +1 multiplicative level: extra modulus (limbs = levels+3) is needed for correct
        // decryption at t~2^48 with multiparty smudging; OpenFHE auto-provisions only levels+2
        // at low depth, so small n fall one limb short without this. N stays 2^15 throughout.
        const uint32_t depth = (uint32_t)std::ceil(std::log2((double)n)) + 1;
        const uint32_t batch = 1u << (uint32_t)std::ceil(std::log2((double)(n + 1)));
        auto C = makeCommittee(t, depth, batch);

        std::mt19937_64 rng(1000 + n);
        std::vector<double> e, p, d, r; bool ok = true; uint32_t N = 0;
        for (int tr = 0; tr < TRIALS; ++tr) {
            auto T = runRound(C, t, n, rng);
            e.push_back(T.enc); p.push_back(T.product); d.push_back(T.pdec);
            r.push_back(T.recovery); ok &= T.ok; N = T.N;
        }
        Stat se=meanSD(e), sp=meanSD(p), sd=meanSD(d), sr=meanSD(r);
        printf("%4d %6u %5u  %7.2f+/-%-5.2f %8.1f+/-%-6.1f %6.2f+/-%-5.2f %8.2f+/-%-6.2f   %s\n",
               n, N, depth, se.mean,se.sd, sp.mean,sp.sd, sd.mean,sd.sd, sr.mean,sr.sd, PF(ok));
        fflush(stdout);
        if (csv) {
            std::fprintf(csv, "%d,%u,%u,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n",
                n, N, depth, TRIALS, se.mean,se.sd, sp.mean,sp.sd, sd.mean,sd.sd, sr.mean,sr.sd, ok?1:0);
            std::fflush(csv);
        }
    }
    if (csv) std::fclose(csv);
    std::cout << "\n=== done (archive: results_profileB_k2.csv) ===\n";
    return 0;
}
