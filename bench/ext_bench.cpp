// Feasibility benchmark: PnF over F_{t^2} (moderate-size messages), swept over n.
// Goal: show the COST of the extension vs. the base single-field (F_t) construction.
// Each F_{t^2} value = two F_t ciphertexts; homomorphic product uses the schoolbook
// F_{t^2} multiply (c = smallest non-residue). Depth = ceil(log2 n) + 3 (leeway, verified).
// Reports F_{t^2} product + recovery (mean +/- SD over TRIALS); recovery = interpolate +
// FindRoots over F_{t^2} (NTL ZZ_pE). Correctness checked every trial (decrypt + multiset).
// This is a feasibility exposition, NOT an optimized implementation (optimization = future work).
#include "openfhe.h"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <iostream>
#include <vector>
#include <map>
#include <random>
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
    const int TRIALS=(argc>1)?std::atoi(argv[1]):3;
    NTL::ZZ Mmod=NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start=((NTL::conv<NTL::ZZ>((long)1)<<48)/Mmod)*Mmod+1;
    while(!NTL::ProbPrime(start)) start+=Mmod;
    const int64_t t=NTL::conv<long>(start);
    auto norm=[&](int64_t v){return((v%t)+t)%t;};

    // F_{t^2}=F_t[Y]/(Y^2-c), c smallest quadratic non-residue (small -> cheap reduction)
    NTL::ZZ_p::init(NTL::conv<NTL::ZZ>((long)t));
    NTL::ZZ pexp=(NTL::conv<NTL::ZZ>((long)t)-1)/2;
    long cval=2; { NTL::ZZ_p cand; while(true){cand=NTL::conv<NTL::ZZ_p>(cval);if(NTL::power(cand,pexp)==NTL::conv<NTL::ZZ_p>((long)-1))break;cval++;} }
    { NTL::ZZ_pX g; NTL::SetCoeff(g,2,1); NTL::SetCoeff(g,0,-NTL::conv<NTL::ZZ_p>(cval)); NTL::ZZ_pE::init(g); }

    printf("=== PnF over F_{t^2} (moderate messages): cost vs n  [t~2^48, c=%ld, %d trials] ===\n",cval,TRIALS);
    printf("   n   depth      N  limbs     product (s)        recovery (ms)     result\n");
    FILE* csv=std::fopen("results_extF2.csv","w");
    if(csv) std::fprintf(csv,"n,depth,N,limbs,product_mean_s,product_sd_s,recovery_mean_ms,recovery_sd_ms,ok\n");

    for(int n : {50,100,200,300,500}){
        const uint32_t depth=(uint32_t)std::ceil(std::log2((double)n))+(n<=200?3u:4u);   // fine-tuned minimal correct depth per n (n<=200: +3/N=2^15; n>=300: +4/N=2^16); verified
        const uint32_t batch=1u<<(uint32_t)std::ceil(std::log2((double)(n+1)));
        CCParams<CryptoContextBFVRNS> p;
        p.SetPlaintextModulus(t); p.SetMultiplicativeDepth(depth);
        p.SetSecurityLevel(HEStd_128_classic); p.SetBatchSize(batch);
        auto cc=GenCryptoContext(p);
        cc->Enable(PKE);cc->Enable(KEYSWITCH);cc->Enable(LEVELEDSHE);
        auto kp=cc->KeyGen(); cc->EvalMultKeyGen(kp.secretKey);
        uint32_t N=cc->GetRingDimension(); size_t limbs=cc->GetCryptoParameters()->GetElementParams()->GetParams().size();
        const long D=n+1;
        std::vector<int64_t> clane(D,cval); auto cpt=cc->MakePackedPlaintext(clane);

        std::vector<double> tp,tr; bool ok=true;
        for(int tr_=0;tr_<TRIALS;++tr_){
            std::mt19937_64 rng(1000*n+tr_);
            std::vector<long> r0(n),r1(n);
            for(int i=0;i<n;i++){r0[i]=(long)(rng()%(uint64_t)t);r1[i]=(long)(rng()%(uint64_t)t);}
            std::vector<Ciphertext<DCRTPoly>> A0(n),A1(n);
            for(int i=0;i<n;i++){std::vector<int64_t> l0(D),l1(D);for(long j=0;j<D;j++){l0[j]=norm((j+1)-r0[i]);l1[j]=norm(-r1[i]);}A0[i]=cc->Encrypt(kp.publicKey,cc->MakePackedPlaintext(l0));A1[i]=cc->Encrypt(kp.publicKey,cc->MakePackedPlaintext(l1));}
            // F_{t^2} product tree (timed)
            auto a=Clock::now();
            std::vector<Ciphertext<DCRTPoly>> P0=A0,P1=A1;
            while(P0.size()>1){std::vector<Ciphertext<DCRTPoly>> Q0,Q1;
                for(size_t i=0;i+1<P0.size();i+=2){auto a0b0=cc->EvalMult(P0[i],P0[i+1]);auto a1b1=cc->EvalMult(P1[i],P1[i+1]);auto a0b1=cc->EvalMult(P0[i],P1[i+1]);auto a1b0=cc->EvalMult(P1[i],P0[i+1]);Q0.push_back(cc->EvalAdd(a0b0,cc->EvalMult(a1b1,cpt)));Q1.push_back(cc->EvalAdd(a0b1,a1b0));}
                if(P0.size()%2==1){Q0.push_back(P0.back());Q1.push_back(P1.back());}P0=Q0;P1=Q1;}
            auto b=Clock::now();
            Plaintext d0,d1;cc->Decrypt(kp.secretKey,P0[0],&d0);cc->Decrypt(kp.secretKey,P1[0],&d1);d0->SetLength(D);d1->SetLength(D);
            auto m0=d0->GetPackedValue(),m1=d1->GetPackedValue();
            // decrypt-correctness guard: expected M(a_j)=prod_i (a_j - r_i) over F_{t^2}.
            // If decryption is garbage we SKIP FindRoots (which would otherwise spin forever on a
            // non-splitting polynomial) and mark the trial FAIL -- this never hangs.
            bool okEval=true;
            for(long j=0;j<D && okEval;j++){
                NTL::ZZ_pE aj=NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_p>((long)(j+1)));
                NTL::ZZ_pE pr; NTL::conv(pr,1L);
                for(int i=0;i<n;i++){NTL::ZZ_pX ri;NTL::SetCoeff(ri,0,NTL::conv<NTL::ZZ_p>((long)r0[i]));NTL::SetCoeff(ri,1,NTL::conv<NTL::ZZ_p>((long)r1[i]));pr*=(aj-NTL::conv<NTL::ZZ_pE>(ri));}
                NTL::ZZ_pX gp;NTL::SetCoeff(gp,0,NTL::conv<NTL::ZZ_p>((long)norm(m0[j])));NTL::SetCoeff(gp,1,NTL::conv<NTL::ZZ_p>((long)norm(m1[j])));
                if(NTL::conv<NTL::ZZ_pE>(gp)!=pr) okEval=false;
            }
            // recovery (timed): interpolate + root-find over F_{t^2}, only if decryption is correct
            auto cstart=Clock::now();
            NTL::vec_ZZ_pE rts;
            if(okEval){
                NTL::vec_ZZ_pE xa,ya;xa.SetLength(D);ya.SetLength(D);
                for(long j=0;j<D;j++){xa[j]=NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_p>((long)(j+1)));NTL::ZZ_pX e;NTL::SetCoeff(e,0,NTL::conv<NTL::ZZ_p>((long)norm(m0[j])));NTL::SetCoeff(e,1,NTL::conv<NTL::ZZ_p>((long)norm(m1[j])));ya[j]=NTL::conv<NTL::ZZ_pE>(e);}
                NTL::ZZ_pEX M;NTL::interpolate(M,xa,ya);NTL::MakeMonic(M);
                NTL::FindRoots(rts,M);
            }
            auto cend=Clock::now();
            // correctness: multiset of (r0,r1) recovered
            std::map<std::pair<long,long>,int> got,exp2;
            for(long i=0;i<rts.length();i++){NTL::ZZ_pX e=NTL::rep(rts[i]);long a0=NTL::deg(e)>=0?NTL::conv<long>(NTL::rep(NTL::coeff(e,0))):0;long a1=NTL::deg(e)>=1?NTL::conv<long>(NTL::rep(NTL::coeff(e,1))):0;got[{a0,a1}]++;}
            for(int i=0;i<n;i++) exp2[{r0[i],r1[i]}]++;
            bool good=okEval&&(got==exp2)&&((long)rts.length()==n);
            ok&=good;
            tp.push_back(ms(a,b)/1000.0); tr.push_back(ms(cstart,cend));
        }
        Stat sp=meanSD(tp),sr=meanSD(tr);
        printf("%4d   %5u  %6u  %5zu   %7.2f+/-%-6.2f   %8.1f+/-%-7.1f  %s\n",n,depth,N,limbs,sp.mean,sp.sd,sr.mean,sr.sd,ok?"PASS":"FAIL");
        if(csv){std::fprintf(csv,"%d,%u,%u,%zu,%.4f,%.4f,%.4f,%.4f,%d\n",n,depth,N,limbs,sp.mean,sp.sd,sr.mean,sr.sd,ok?1:0);std::fflush(csv);}
    }
    if(csv)std::fclose(csv);
    printf("=== done (archive: results_extF2.csv) ===\n");
    return 0;
}
