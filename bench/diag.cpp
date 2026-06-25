// Per-phase timing diagnostic for a single n=50 round at t~2^48, to locate the regression.
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
using namespace lbcrypto;
using Clock = std::chrono::steady_clock;
static double ms(Clock::time_point a, Clock::time_point b){ return std::chrono::duration<double,std::milli>(b-a).count(); }
#define MARK(label, ...) { auto _a=Clock::now(); __VA_ARGS__; auto _b=Clock::now(); printf("%-18s %10.1f ms\n", label, ms(_a,_b)); fflush(stdout); }

int main(int argc, char** argv){
    int n = (argc>1)?std::atoi(argv[1]):50;
    int tbits = (argc>2)?std::atoi(argv[2]):48;
    int dpad  = (argc>3)?std::atoi(argv[3]):0;
    NTL::ZZ Mmod = NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start = ((NTL::conv<NTL::ZZ>((long)1) << tbits) / Mmod) * Mmod + 1;
    while(!NTL::ProbPrime(start)) start += Mmod;
    const int64_t t = NTL::conv<long>(start);
    const uint32_t depth = (uint32_t)std::ceil(std::log2((double)n)) + (uint32_t)dpad;
    const uint32_t batch = 1u << (uint32_t)std::ceil(std::log2((double)(n+1)));
    printf("n=%d t=%lld (~2^48) depth=%u\n", n, (long long)t, depth); fflush(stdout);

    CryptoContext<DCRTPoly> cc;
    MARK("context", {
        CCParams<CryptoContextBFVRNS> p;
        p.SetPlaintextModulus(t); p.SetMultiplicativeDepth(depth);
        p.SetSecurityLevel(HEStd_128_classic); p.SetBatchSize(batch);
        cc = GenCryptoContext(p);
        cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(MULTIPARTY);
    });
    printf("N=%u limbs=%zu\n", cc->GetRingDimension(),
           cc->GetCryptoParameters()->GetElementParams()->GetParams().size()); fflush(stdout);

    KeyPair<DCRTPoly> kp1, kp2;
    MARK("keygen", {
        kp1 = cc->KeyGen(); kp2 = cc->MultipartyKeyGen(kp1.publicKey);
        auto emk=cc->KeySwitchGen(kp1.secretKey,kp1.secretKey);
        auto emk2=cc->MultiKeySwitchGen(kp2.secretKey,kp2.secretKey,emk);
        auto emkAB=cc->MultiAddEvalKeys(emk,emk2,kp2.publicKey->GetKeyTag());
        auto emkBAB=cc->MultiMultEvalKey(kp2.secretKey,emkAB,kp2.publicKey->GetKeyTag());
        auto emkAAB=cc->MultiMultEvalKey(kp1.secretKey,emkAB,kp2.publicKey->GetKeyTag());
        auto emkF=cc->MultiAddEvalMultKeys(emkAAB,emkBAB,emkAB->GetKeyTag());
        cc->InsertEvalMultKey({emkF});
    });
    auto jointPk = kp2.publicKey;
    const long D = n+1;
    auto norm=[&](int64_t v){ return ((v%t)+t)%t; };
    std::mt19937_64 rng(1050);
    std::vector<int64_t> roots(n); std::map<long,int> seen;
    for(int i=0;i<n;++i){ int64_t r; do{ r=(int64_t)(rng()%(uint64_t)t);}while(seen[(long)r]++); roots[i]=r; }

    std::vector<Ciphertext<DCRTPoly>> cts;
    MARK("encrypt(n)", {
        for(int i=0;i<n;++i){ std::vector<int64_t> lanes(D); for(long j=0;j<D;++j) lanes[j]=norm((j+1)-roots[i]);
            cts.push_back(cc->Encrypt(jointPk, cc->MakePackedPlaintext(lanes))); }
    });
    Ciphertext<DCRTPoly> prod;
    MARK("product", {
        while(cts.size()>1){ std::vector<Ciphertext<DCRTPoly>> nx;
            for(size_t i=0;i+1<cts.size();i+=2) nx.push_back(cc->EvalMult(cts[i],cts[i+1]));
            if(cts.size()%2==1) nx.push_back(cts.back()); cts=nx; }
        prod=cts[0];
    });
    MARK("sanitize", { std::vector<int64_t> zl(D,0); prod=cc->EvalAdd(prod, cc->Encrypt(jointPk, cc->MakePackedPlaintext(zl))); });
    std::vector<Ciphertext<DCRTPoly>> pL,pM;
    MARK("pdec lead", { pL=cc->MultipartyDecryptLead({prod}, kp1.secretKey); });
    MARK("pdec main", { pM=cc->MultipartyDecryptMain({prod}, kp2.secretKey); });
    Plaintext res; std::vector<int64_t> evals;
    MARK("combine", { cc->MultipartyDecryptFusion({pL[0],pM[0]}, &res); res->SetLength(D); evals=res->GetPackedValue(); for(auto&v:evals)v=norm(v); });
    bool okEval=true; for(long j=0;j<D && okEval;++j){ __int128 e=1; for(int i=0;i<n;++i) e=e*norm((j+1)-roots[i])%t; if((int64_t)e!=evals[j]) okEval=false; }
    printf("decrypt correct:   %s\n", okEval?"YES":"NO"); fflush(stdout);
    if(!okEval){ printf("(skipping findroots: M is garbage)\n"); return 2; }
    NTL::ZZ_pX M; NTL::vec_ZZ_p rts;
    MARK("interpolate", {
        NTL::ZZ_p::init(NTL::conv<NTL::ZZ>((long)t));
        NTL::vec_ZZ_p xa,ya; xa.SetLength(D); ya.SetLength(D);
        for(long j=0;j<D;++j){ xa[j]=NTL::conv<NTL::ZZ_p>((long)(j+1)); ya[j]=NTL::conv<NTL::ZZ_p>((long)evals[j]); }
        NTL::interpolate(M,xa,ya); NTL::MakeMonic(M);
    });
    MARK("findroots", { NTL::FindRoots(rts,M); });
    printf("roots found: %ld (expected %d)\n", rts.length(), n); fflush(stdout);
    return 0;
}
