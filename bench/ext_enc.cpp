// Per-sender encryption cost for the F_{t^2} (longer-message) construction: a sender encrypts
// TWO ordinary F_t ciphertexts (the two coordinates of its F_{t^2} value). Reports enc/sender
// (= time for the 2 encryptions) per n, mean +/- SD. Same context as ext_bench (depth ceil(log2 n)+3).
#include "openfhe.h"
#include <NTL/ZZ.h>
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
using namespace lbcrypto;
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
struct Stat{double mean,sd;};
static Stat meanSD(const std::vector<double>&v){double m=0;for(double x:v)m+=x;m/=v.size();double s=0;for(double x:v)s+=(x-m)*(x-m);s=v.size()>1?std::sqrt(s/(v.size()-1)):0;return{m,s};}

int main(int argc,char**argv){
    setvbuf(stdout,NULL,_IONBF,0);
    const int TRIALS=(argc>1)?std::atoi(argv[1]):10;
    NTL::ZZ Mmod=NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start=((NTL::conv<NTL::ZZ>((long)1)<<48)/Mmod)*Mmod+1;
    while(!NTL::ProbPrime(start)) start+=Mmod;
    const int64_t t=NTL::conv<long>(start);
    auto norm=[&](int64_t v){return((v%t)+t)%t;};
    printf("=== F_{t^2} per-sender encryption (2 ciphertexts), %d trials ===\n",TRIALS);
    printf("   n  limbs   enc/sender (ms)\n");
    for(int n : {50,100,200,300,500}){
        const uint32_t depth=(uint32_t)std::ceil(std::log2((double)n))+(n<=200?3u:4u);   // match ext_bench (fine-tuned per-n)
        const uint32_t batch=1u<<(uint32_t)std::ceil(std::log2((double)(n+1)));
        CCParams<CryptoContextBFVRNS> p;
        p.SetPlaintextModulus(t); p.SetMultiplicativeDepth(depth);
        p.SetSecurityLevel(HEStd_128_classic); p.SetBatchSize(batch);
        auto cc=GenCryptoContext(p); cc->Enable(PKE);cc->Enable(KEYSWITCH);cc->Enable(LEVELEDSHE);
        auto kp=cc->KeyGen();
        size_t limbs=cc->GetCryptoParameters()->GetElementParams()->GetParams().size();
        const long D=n+1;
        std::vector<int64_t> l0(D),l1(D); for(long j=0;j<D;j++){l0[j]=norm((j+1)-12345);l1[j]=norm(-6789);}
        auto pt0=cc->MakePackedPlaintext(l0), pt1=cc->MakePackedPlaintext(l1);
        std::vector<double> e;
        for(int tr=0;tr<TRIALS;++tr){ auto a=Clock::now(); auto c0=cc->Encrypt(kp.publicKey,pt0); auto c1=cc->Encrypt(kp.publicKey,pt1); auto b=Clock::now(); (void)c0;(void)c1; e.push_back(ms(a,b)); }
        Stat se=meanSD(e);
        printf("%4d  %5zu   %6.2f+/-%-5.2f\n",n,limbs,se.mean,se.sd); fflush(stdout);
    }
    return 0;
}
