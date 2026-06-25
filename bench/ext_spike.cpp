// Feasibility spike: PnF over an extension field F_{t^s} (here s=2) to enlarge the message domain.
// Each F_{t^2} value is carried by TWO ordinary F_t ciphertexts (its coordinates over F_t); the
// homomorphic product uses the schoolbook F_{t^2} multiply (a0+a1 Y)(b0+b1 Y) = (a0b0 + c a1b1) +
// (a0b1 + a1b0) Y, with F_{t^2}=F_t[Y]/(Y^2 - c) and c the SMALLEST quadratic non-residue (tiny, so
// its scalar multiply adds negligible noise -> same FHE params as the base, ~s^2x the product cost).
// Recovery: interpolate + root-find over F_{t^2} via NTL ZZ_pE. Verifies multiset + tag recovery.
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
#include <cmath>
#include <cstdint>
#include <chrono>
using namespace lbcrypto;
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
static const char* PF(bool ok){ return ok?"PASS":"FAIL"; }

int main(int argc, char** argv){
    setvbuf(stdout,NULL,_IONBF,0);
    int n = (argc>1)? std::atoi(argv[1]) : 8;
    int dpad = (argc>2)? std::atoi(argv[2]) : 0;
    NTL::ZZ Mmod = NTL::conv<NTL::ZZ>((long)131072);
    NTL::ZZ start = ((NTL::conv<NTL::ZZ>((long)1)<<48)/Mmod)*Mmod + 1;
    while(!NTL::ProbPrime(start)) start += Mmod;
    const int64_t t = NTL::conv<long>(start);
    const uint32_t depth = (uint32_t)std::ceil(std::log2((double)n)) + 1 + (uint32_t)dpad;
    const uint32_t batch = 1u << (uint32_t)std::ceil(std::log2((double)(n+1)));

    // --- NTL F_{t^2} = F_t[Y]/(Y^2 - c), c = smallest quadratic non-residue ---
    NTL::ZZ_p::init(NTL::conv<NTL::ZZ>((long)t));
    NTL::ZZ exp = (NTL::conv<NTL::ZZ>((long)t)-1)/2;
    long cval=2; NTL::ZZ_p c;
    while(true){ NTL::ZZ_p cand=NTL::conv<NTL::ZZ_p>(cval); if(NTL::power(cand,exp)==NTL::conv<NTL::ZZ_p>((long)-1)){c=cand;break;} cval++; }
    NTL::ZZ_pX g; NTL::SetCoeff(g,2,1); NTL::SetCoeff(g,0,-c);
    NTL::ZZ_pE::init(g);
    printf("ext_spike: F_{t^2}, t=%lld (~2^48), c=%ld (smallest non-residue), n=%d, depth=%u\n",(long long)t,cval,n,depth);

    // --- OpenFHE BFV (standard F_t packed) ---
    auto ta=Clock::now();
    CCParams<CryptoContextBFVRNS> p;
    p.SetPlaintextModulus(t); p.SetMultiplicativeDepth(depth);
    p.SetSecurityLevel(HEStd_128_classic); p.SetBatchSize(batch);
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
    auto kp = cc->KeyGen(); cc->EvalMultKeyGen(kp.secretKey);
    printf("N=%u, limbs=%zu\n", cc->GetRingDimension(), cc->GetCryptoParameters()->GetElementParams()->GetParams().size());
    printf("setup       %.0f ms\n", ms(ta,Clock::now())); auto tb=Clock::now();
    auto norm=[&](int64_t v){ return ((v%t)+t)%t; };
    const long D = n+1;

    // --- senders: roots r_i = (r0,r1) in F_{t^2} (each coord a full ~48-bit field element) ---
    std::mt19937_64 rng(2024);
    std::vector<long> r0(n), r1(n);
    for(int i=0;i<n;i++){ r0[i]=(long)(rng()%(uint64_t)t); r1[i]=(long)(rng()%(uint64_t)t); }

    // factor i evaluated at base-field points a_j=j+1: (a_j - r0_i, -r1_i) in F_{t^2}
    auto enc=[&](const std::vector<int64_t>& l){ return cc->Encrypt(kp.publicKey, cc->MakePackedPlaintext(l)); };
    std::vector<Ciphertext<DCRTPoly>> A0(n), A1(n);
    for(int i=0;i<n;i++){
        std::vector<int64_t> l0(D), l1(D);
        for(long j=0;j<D;j++){ l0[j]=norm((j+1)-r0[i]); l1[j]=norm(-r1[i]); }
        A0[i]=enc(l0); A1[i]=enc(l1);
    }
    printf("encrypt(2n)  %.0f ms\n", ms(tb,Clock::now())); auto tc=Clock::now();
    std::vector<int64_t> clane(D, cval); auto cpt = cc->MakePackedPlaintext(clane);
    auto fmul=[&](const Ciphertext<DCRTPoly>&a0,const Ciphertext<DCRTPoly>&a1,
                  const Ciphertext<DCRTPoly>&b0,const Ciphertext<DCRTPoly>&b1,
                  Ciphertext<DCRTPoly>&o0, Ciphertext<DCRTPoly>&o1){
        auto a0b0=cc->EvalMult(a0,b0); auto a1b1=cc->EvalMult(a1,b1);
        auto a0b1=cc->EvalMult(a0,b1); auto a1b0=cc->EvalMult(a1,b0);
        o0=cc->EvalAdd(a0b0, cc->EvalMult(a1b1,cpt));   // a0b0 + c*a1b1
        o1=cc->EvalAdd(a0b1, a1b0);                      // a0b1 + a1b0
    };
    // product tree over F_{t^2}
    std::vector<Ciphertext<DCRTPoly>> P0=A0, P1=A1;
    while(P0.size()>1){
        std::vector<Ciphertext<DCRTPoly>> Q0,Q1;
        for(size_t i=0;i+1<P0.size();i+=2){ Ciphertext<DCRTPoly> o0,o1; fmul(P0[i],P1[i],P0[i+1],P1[i+1],o0,o1); Q0.push_back(o0); Q1.push_back(o1); }
        if(P0.size()%2==1){ Q0.push_back(P0.back()); Q1.push_back(P1.back()); }
        P0=Q0; P1=Q1;
    }
    printf("product      %.0f ms\n", ms(tc,Clock::now())); auto td=Clock::now();
    Plaintext d0,d1; cc->Decrypt(kp.secretKey,P0[0],&d0); cc->Decrypt(kp.secretKey,P1[0],&d1);
    d0->SetLength(D); d1->SetLength(D);
    auto m0=d0->GetPackedValue(), m1=d1->GetPackedValue();
    printf("decrypt      %.0f ms\n", ms(td,Clock::now()));

    // --- decrypt-correctness guard: expected M(a_j)=prod_i (a_j - r_i) over F_{t^2} ---
    bool okEval=true;
    for(long j=0;j<D && okEval;j++){
        NTL::ZZ_pE aj=NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_p>((long)(j+1)));
        NTL::ZZ_pE prod; NTL::conv(prod,1L);
        for(int i=0;i<n;i++){ NTL::ZZ_pX ri; NTL::SetCoeff(ri,0,NTL::conv<NTL::ZZ_p>((long)r0[i])); NTL::SetCoeff(ri,1,NTL::conv<NTL::ZZ_p>((long)r1[i])); prod*=(aj-NTL::conv<NTL::ZZ_pE>(ri)); }
        NTL::ZZ_pX gp; NTL::SetCoeff(gp,0,NTL::conv<NTL::ZZ_p>((long)norm(m0[j]))); NTL::SetCoeff(gp,1,NTL::conv<NTL::ZZ_p>((long)norm(m1[j])));
        if(NTL::conv<NTL::ZZ_pE>(gp)!=prod) okEval=false;
    }
    printf("decrypt correct: %s\n", PF(okEval));
    if(!okEval){ printf("(M garbage -> skip FindRoots; raise dpad)\n"); return 2; }
    auto te=Clock::now();

    // --- recover over F_{t^2}: interpolate M(x), find roots ---
    NTL::vec_ZZ_pE xa, ya; xa.SetLength(D); ya.SetLength(D);
    for(long j=0;j<D;j++){
        xa[j]=NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_p>((long)(j+1)));
        NTL::ZZ_pX e; NTL::SetCoeff(e,0,NTL::conv<NTL::ZZ_p>((long)norm(m0[j]))); NTL::SetCoeff(e,1,NTL::conv<NTL::ZZ_p>((long)norm(m1[j])));
        ya[j]=NTL::conv<NTL::ZZ_pE>(e);
    }
    NTL::ZZ_pEX M; NTL::interpolate(M,xa,ya); NTL::MakeMonic(M);
    NTL::vec_ZZ_pE rts; NTL::FindRoots(rts,M);
    printf("recovery     %.0f ms\n", ms(te,Clock::now()));

    // --- compare multiset of (r0,r1) pairs ---
    auto key=[&](long a,long b){ return std::make_pair(a,b); };
    std::map<std::pair<long,long>,int> got, exp2;
    for(long i=0;i<rts.length();i++){ NTL::ZZ_pX e=NTL::rep(rts[i]); long a=NTL::deg(e)>=0?NTL::conv<long>(NTL::rep(NTL::coeff(e,0))):0; long b=NTL::deg(e)>=1?NTL::conv<long>(NTL::rep(NTL::coeff(e,1))):0; got[key(a,b)]++; }
    for(int i=0;i<n;i++) exp2[key(r0[i],r1[i])]++;
    bool okMs = (got==exp2) && ((long)rts.length()==n);
    // sender verifiability: each sender finds its own (r0,r1)
    bool okSV=true; for(int i=0;i<n;i++) if(!got.count(key(r0[i],r1[i]))) okSV=false;

    printf("[1] roots found: %ld (expected %d)\n", rts.length(), n);
    printf("[2] multiset recovered over F_{t^2}: %s\n", PF(okMs));
    printf("[3] sender-verifiability (each tag/value found): %s\n", PF(okSV));
    printf("=== ext-field feasibility: %s ===\n", PF(okMs&&okSV));
    return (okMs&&okSV)?0:1;
}
