// Probe: largest packed-encoding plaintext modulus t that BFVRNS accepts at depth 10
// (the n=1000 setting), 128-bit security, 64-bit native backend. For each candidate
// 2^b we pick the smallest prime >= 2^b with t = 1 (mod 2^17), build the context,
// and try a pack+encrypt+mult+threshold-decrypt round; report PASS/FAIL, N, and the
// recovered-correctly flag. This tells us the real tag-length ceiling.
#include "openfhe.h"
#include <NTL/ZZ.h>
#include <iostream>
#include <vector>
#include <cstdint>
using namespace lbcrypto;

static int64_t primeAt(int bits) {
    NTL::ZZ Mmod  = NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start = ((NTL::conv<NTL::ZZ>((long)1) << bits) / Mmod) * Mmod + 1;
    while (!NTL::ProbPrime(start)) start += Mmod;
    return NTL::conv<long>(start);
}

int main() {
    const uint32_t depth = 10, batch = 1024;
    std::cout << "bits  t                      build  N      encode/mult/decrypt\n";
    for (int bits : {44, 48, 50, 52, 54, 55, 56, 57, 58, 59, 60}) {
        int64_t t = primeAt(bits);
        try {
            CCParams<CryptoContextBFVRNS> p;
            p.SetPlaintextModulus(t);
            p.SetMultiplicativeDepth(depth);
            p.SetSecurityLevel(HEStd_128_classic);
            p.SetBatchSize(batch);
            auto cc = GenCryptoContext(p);
            cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(MULTIPARTY);
            uint32_t N = cc->GetRingDimension();
            // quick functional check: two-party threshold, product of 4 factors
            auto kp1 = cc->KeyGen();
            auto kp2 = cc->MultipartyKeyGen(kp1.publicKey);
            auto emk   = cc->KeySwitchGen(kp1.secretKey, kp1.secretKey);
            auto emk2  = cc->MultiKeySwitchGen(kp2.secretKey, kp2.secretKey, emk);
            auto emkAB = cc->MultiAddEvalKeys(emk, emk2, kp2.publicKey->GetKeyTag());
            auto emkBAB= cc->MultiMultEvalKey(kp2.secretKey, emkAB, kp2.publicKey->GetKeyTag());
            auto emkAAB= cc->MultiMultEvalKey(kp1.secretKey, emkAB, kp2.publicKey->GetKeyTag());
            auto emkF  = cc->MultiAddEvalMultKeys(emkAAB, emkBAB, emkAB->GetKeyTag());
            cc->InsertEvalMultKey({emkF});
            auto norm = [&](int64_t v){ return ((v % t) + t) % t; };
            const long D = 5;
            std::vector<int64_t> roots = {123456789LL, 987654321LL, 555555555LL, 111111111LL};
            std::vector<Ciphertext<DCRTPoly>> cts;
            for (int i = 0; i < 4; ++i) {
                std::vector<int64_t> lanes(D);
                for (long j = 0; j < D; ++j) lanes[j] = norm((j+1) - roots[i]);
                cts.push_back(cc->Encrypt(kp2.publicKey, cc->MakePackedPlaintext(lanes)));
            }
            auto prod = cc->EvalMult(cc->EvalMult(cts[0],cts[1]), cc->EvalMult(cts[2],cts[3]));
            auto pL = cc->MultipartyDecryptLead({prod}, kp1.secretKey);
            auto pM = cc->MultipartyDecryptMain({prod}, kp2.secretKey);
            Plaintext res; cc->MultipartyDecryptFusion({pL[0],pM[0]}, &res);
            res->SetLength(D); auto ev = res->GetPackedValue();
            bool ok = true;
            for (long j = 0; j < D; ++j) {
                __int128 e = 1;
                for (int i = 0; i < 4; ++i) e = e * norm((j+1)-roots[i]) % t;
                if ((int64_t)e != norm(ev[j])) ok = false;
            }
            printf("%4d  %-22lld build  N=%-5u  %s\n", bits, (long long)t, N, ok?"PASS":"WRONG");
        } catch (const std::exception& e) {
            printf("%4d  %-22lld FAIL   (%s)\n", bits, (long long)t, e.what());
        }
        fflush(stdout);
    }
    return 0;
}
