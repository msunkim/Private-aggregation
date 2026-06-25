// Parallel-scaling of the aggregator's homomorphic product vs OpenMP thread count.
// Fixes n (default 300), builds the product inputs once, then times the product tree at
// thread counts {1,2,4,8,12} (omp_set_num_threads between runs), reporting mean+/-SD and
// speedup over single-threaded. The product is the only n-growing cost; this substantiates
// that it parallelizes on the cloud. (Single-key context: EvalMult cost is key-count
// independent, so this matches the threshold setting's product.)
#include "openfhe.h"
extern "C" void omp_set_num_threads(int);   // from libomp (linked); avoids needing omp.h on the include path
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <NTL/ZZ.h>
using namespace lbcrypto;
using Clock = std::chrono::steady_clock;
static double ms(Clock::time_point a, Clock::time_point b){ return std::chrono::duration<double,std::milli>(b-a).count(); }

int main(int argc, char** argv){
    const int n = (argc>1)?std::atoi(argv[1]):300;
    const int TRIALS = (argc>2)?std::atoi(argv[2]):3;
    const int only = (argc>3)?std::atoi(argv[3]):0;   // if >0, measure ONLY this thread count (fresh process)
    NTL::ZZ Mmod = NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start = ((NTL::conv<NTL::ZZ>((long)1) << 48) / Mmod) * Mmod + 1;
    while(!NTL::ProbPrime(start)) start += Mmod;
    const int64_t t = NTL::conv<long>(start);
    const uint32_t depth = (uint32_t)std::ceil(std::log2((double)n)) + 1;
    const uint32_t batch = 1u << (uint32_t)std::ceil(std::log2((double)(n+1)));

    CCParams<CryptoContextBFVRNS> p;
    p.SetPlaintextModulus(t); p.SetMultiplicativeDepth(depth);
    p.SetSecurityLevel(HEStd_128_classic); p.SetBatchSize(batch);
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
    auto kp = cc->KeyGen(); cc->EvalMultKeyGen(kp.secretKey);
    auto norm=[&](int64_t v){ return ((v%t)+t)%t; };
    const long D = n+1;

    std::vector<Ciphertext<DCRTPoly>> base;
    for(int i=0;i<n;++i){ std::vector<int64_t> lanes(D); for(long j=0;j<D;++j) lanes[j]=norm((j+1)-(1000+i));
        base.push_back(cc->Encrypt(kp.publicKey, cc->MakePackedPlaintext(lanes))); }

    auto product=[&](std::vector<Ciphertext<DCRTPoly>> cts){
        while(cts.size()>1){ std::vector<Ciphertext<DCRTPoly>> nx;
            for(size_t i=0;i+1<cts.size();i+=2) nx.push_back(cc->EvalMult(cts[i],cts[i+1]));
            if(cts.size()%2==1) nx.push_back(cts.back()); cts=nx; } };

    printf("=== Product parallel-scaling (n=%d, depth=%u, N=%u, t~2^48, %d trials) ===\n",
           n, depth, cc->GetRingDimension(), TRIALS);
    printf("threads   product (s)        speedup\n");
    double base1=0;
    std::vector<int> tcs = (only>0) ? std::vector<int>{only} : std::vector<int>{1,2,4,8,12};
    for(int tc : tcs){
        omp_set_num_threads(tc);
        std::vector<double> v;
        for(int tr=0;tr<TRIALS;++tr){ auto a=Clock::now(); product(base); auto b=Clock::now(); v.push_back(ms(a,b)/1000.0); }
        double m=0; for(double x:v)m+=x; m/=v.size(); double s=0; for(double x:v)s+=(x-m)*(x-m); s=v.size()>1?std::sqrt(s/(v.size()-1)):0;
        if(tc==1) base1=m;
        printf("%5d    %7.2f+/-%-5.2f    %5.2fx\n", tc, m, s, base1/m);
        fflush(stdout);
    }
    return 0;
}
