// Ciphertext / transmission-size probe for the Communication paragraph.
// For each n it builds the same BFV context as profile_b_bench, reads the RNS limb count,
// and reports the byte sizes of the objects each party transmits:
//   sender    : one fresh ciphertext (2 ring elements at full modulus)
//   aggregator: one aggregate ciphertext (2 ring elements)
//   committee : one partial decryption (1 ring element)
// (Proof sizes are reported separately by sigma_bench.)

#include "openfhe.h"
#include <NTL/ZZ.h>
#include <iostream>
#include <cmath>
#include <cstdint>
using namespace lbcrypto;

int main() {
    NTL::ZZ Mmod = NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start = ((NTL::conv<NTL::ZZ>((long)1) << 48) / Mmod) * Mmod + 1;
    while (!NTL::ProbPrime(start)) start += Mmod;
    const int64_t t = NTL::conv<long>(start);

    printf("   n  depth      N  limbs  log2q   ciphertext   partial-dec\n");
    for (int n : {50,100,200,300,500,1000}) {
        uint32_t depth = (uint32_t)std::ceil(std::log2((double)n)) + 1;  // +1 headroom (matches profile_b @ t~2^48)
        uint32_t batch = 1u << (uint32_t)std::ceil(std::log2((double)(n+1)));
        CCParams<CryptoContextBFVRNS> p;
        p.SetPlaintextModulus(t);
        p.SetMultiplicativeDepth(depth);
        p.SetSecurityLevel(HEStd_128_classic);
        p.SetBatchSize(batch);
        auto cc = GenCryptoContext(p);
        auto ep = cc->GetCryptoParameters()->GetElementParams();
        size_t limbs = ep->GetParams().size();
        uint32_t N = cc->GetRingDimension();
        uint32_t logq = ep->GetModulus().GetMSB();          // bit-length of the composite ciphertext modulus Q
        double poly = (double)limbs * N * 8.0;            // one ring element, 8 bytes / RNS coeff
        double ct   = 2.0 * poly;                          // ciphertext = 2 ring elements
        printf("%4d  %5u  %6u  %5zu  %5u   %8.2f MB  %8.2f MB\n",
               n, depth, N, limbs, logq, ct/1e6, poly/1e6);
    }
    return 0;
}
