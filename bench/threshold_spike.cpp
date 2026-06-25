// Threshold BFV spike (Profile B / evaluation domain).
// ---------------------------------------------------------------------------------------
// Purpose: validate OpenFHE's NATIVE multiparty (threshold) BFV before building the full
// §6 benchmark. We exercise the complete threshold path that PnF needs:
//
//   1. Multiparty key generation  -> a JOINT public key (no single party holds the secret).
//   2. Multiparty relinearization -> joint eval-mult key, so EvalMult works under the joint key.
//   3. Homomorphic product of n linear factors encrypted under the joint pk (Profile B lanes).
//   4. THRESHOLD decryption: MultipartyDecryptLead + MultipartyDecryptMain + Fusion.
//   5. Interpolate + factor over F_t -> recover the multiset, and check it matches.
//
// Committee size k is fixed to 2 here (lead + one main), the canonical OpenFHE multiparty
// pattern; the eval-mult-key chaining generalizes to larger k. If this passes, threshold
// BFV is sound for our setting and we can scale it up (k, n) in the benchmark proper.
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
#include <cmath>
#include <cstdint>

using namespace lbcrypto;
static const char* PF(bool ok) { return ok ? "PASS" : "FAIL"; }

int main() {
    const int n = 4;                                                   // factors in the product
    const uint32_t depth = (uint32_t)std::ceil(std::log2((double)n));  // multiplicative depth

    // Profile B: large prime t ~ 2^48 with t = 1 (mod 2^17) so packing works for N <= 65536.
    NTL::ZZ Mmod  = NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start = ((NTL::conv<NTL::ZZ>((long)1) << 48) / Mmod) * Mmod + 1;
    while (!NTL::ProbPrime(start)) start += Mmod;
    const int64_t t = NTL::conv<long>(start);

    CCParams<CryptoContextBFVRNS> p;
    p.SetPlaintextModulus(t);
    p.SetMultiplicativeDepth(depth);
    p.SetSecurityLevel(HEStd_128_classic);
    p.SetBatchSize(8);
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(MULTIPARTY);

    std::cout << "=== Threshold BFV spike (k=2, Profile B) ===\n";
    std::cout << "t = " << t << " (~2^48), N = " << cc->GetRingDimension()
              << ", depth = " << depth << ", n = " << n << "\n\n";

    // ---- 1. Multiparty key generation: joint public key ----
    auto kp1 = cc->KeyGen();
    auto kp2 = cc->MultipartyKeyGen(kp1.publicKey);   // chains onto party 1
    auto jointPk = kp2.publicKey;                     // the JOINT public key

    // ---- 2. Multiparty relinearization (joint eval-mult key) ----
    auto evalMultKey   = cc->KeySwitchGen(kp1.secretKey, kp1.secretKey);
    auto evalMultKey2  = cc->MultiKeySwitchGen(kp2.secretKey, kp2.secretKey, evalMultKey);
    auto evalMultAB    = cc->MultiAddEvalKeys(evalMultKey, evalMultKey2, kp2.publicKey->GetKeyTag());
    auto evalMultBAB   = cc->MultiMultEvalKey(kp2.secretKey, evalMultAB, kp2.publicKey->GetKeyTag());
    auto evalMultAAB   = cc->MultiMultEvalKey(kp1.secretKey, evalMultAB, kp2.publicKey->GetKeyTag());
    auto evalMultFinal = cc->MultiAddEvalMultKeys(evalMultAAB, evalMultBAB, evalMultAB->GetKeyTag());
    cc->InsertEvalMultKey({evalMultFinal});

    // ---- 3. Encode n factors (x - r_i) as Profile-B lanes, encrypt under the JOINT pk ----
    const long degM = n;
    const long D    = degM + 1;                       // lanes = evaluation points a_j = j+1
    auto norm = [&](int64_t v) { return ((v % t) + t) % t; };

    std::mt19937_64 rng(7);
    std::vector<int64_t> roots(n);
    for (int i = 0; i < n; ++i) roots[i] = (int64_t)(rng() % (1u << 20));

    std::vector<Ciphertext<DCRTPoly>> cts;
    for (int i = 0; i < n; ++i) {
        std::vector<int64_t> lanes(D);
        for (long j = 0; j < D; ++j) lanes[j] = norm((j + 1) - roots[i]);
        Plaintext pt = cc->MakePackedPlaintext(lanes);
        cts.push_back(cc->Encrypt(jointPk, pt));      // encrypt under the joint public key
    }

    // homomorphic product tree (uses the joint eval-mult key inserted above)
    while (cts.size() > 1) {
        std::vector<Ciphertext<DCRTPoly>> nxt;
        for (size_t i = 0; i + 1 < cts.size(); i += 2) nxt.push_back(cc->EvalMult(cts[i], cts[i + 1]));
        if (cts.size() % 2 == 1) nxt.push_back(cts.back());
        cts = nxt;
    }
    auto prod = cts[0];

    // ---- 4. Threshold (multiparty) decryption ----
    auto partial1 = cc->MultipartyDecryptLead({prod}, kp1.secretKey);
    auto partial2 = cc->MultipartyDecryptMain({prod}, kp2.secretKey);
    Plaintext res;
    cc->MultipartyDecryptFusion({partial1[0], partial2[0]}, &res);
    res->SetLength(D);
    auto evals = res->GetPackedValue();
    for (auto& v : evals) v = norm(v);

    // expected M(a_j) = prod_i (a_j - r_i)
    bool okEval = true;
    for (long j = 0; j < D; ++j) {
        int64_t e = 1;
        for (int i = 0; i < n; ++i) e = (int64_t)((__int128)e * norm((j + 1) - roots[i]) % t);
        if (e != evals[j]) okEval = false;
    }
    std::cout << "[1] joint-key product == expected over " << D << " lanes: " << PF(okEval) << "\n";

    // ---- 5. Interpolate + factor over F_t -> recover the multiset ----
    NTL::ZZ_p::init(NTL::conv<NTL::ZZ>((long)t));
    NTL::vec_ZZ_p xa, ya; xa.SetLength(D); ya.SetLength(D);
    for (long j = 0; j < D; ++j) {
        xa[j] = NTL::conv<NTL::ZZ_p>((long)(j + 1));
        ya[j] = NTL::conv<NTL::ZZ_p>((long)evals[j]);
    }
    NTL::ZZ_pX M; NTL::interpolate(M, xa, ya);
    NTL::MakeMonic(M);
    NTL::vec_ZZ_p rts; NTL::FindRoots(rts, M);
    std::map<long, int> got, exp;
    for (long j = 0; j < rts.length(); ++j) got[NTL::conv<long>(NTL::rep(rts[j]))]++;
    for (int i = 0; i < n; ++i) exp[(long)roots[i]]++;
    bool okMs = (got == exp);
    std::cout << "[2] threshold-decrypt + factor recovers multiset: " << PF(okMs) << "\n";

    std::cout << "\n=== thresholdness: " << PF(okEval && okMs) << " ===\n";
    return (okEval && okMs) ? 0 : 1;
}
