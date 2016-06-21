#define SCHEME "2"
#define NUMMASSES "3"

*
* load default values for preprocessor variables
*

*
* global definition file for boundary values used in the computation
*

*
* default value for DALA1 and TOPOLOGY
*

#ifndef `DALA1'
        #define DALA1 "1"
#endif

#ifndef `TOPOLOGY'
        #define TOPOLOGY "arb"
#endif

*
* default value for the number of LOOPS
*

#ifndef `LOOPS'
        #define LOOPS "0"
#endif

*
* INTLOOPS is used in the declaration and has a minimum value of 1
*

#if `LOOPS' == "0" 
        #define INTLOOPS "1"
        #else
        #define INTLOOPS "`LOOPS'"
#endif

*
* MAXLOOPS is used in cutep 
*

#ifndef `MAXLOOPS'
        #define MAXLOOPS "3"
#endif

*
* default value for the number of external momenta
*

#ifndef `NUMEXTMOM'
        #define NUMEXTMOM "1"
#endif 

*
* default value for the number of masses
*

#ifndef `NUMMASSES'
        #define NUMMASSES "3"
#endif

*
* default value of variants for the indices mu, nu, al, be, is
* 

#ifndef `NUMINDEX'
        #define NUMINDEX "40"
#endif

*
* default depth of pochhammer table
*

#ifndef `POCHHAMMER'
        #define POCHHAMMER "24"
#endif

*
* default maximal number of fermion traces
*

#ifndef `NUMFERMIONTRACES'
        #define NUMFERMIONTRACES "4"
#endif

*
* default SCHEME for MINCER
*

#ifndef `SCHEME'
        #define SCHEME "2"
#endif


*
* global declare file
*


S n,ep(:6);
dimension n;

S eQ,eQ1,...,eQ`NUMEXTMOM';
V q,q1,...,q`NUMEXTMOM', Q,Q1,...,Q`NUMEXTMOM';

S M,m,M1,...,M`NUMMASSES', eM1,...,eM`NUMMASSES';

V P,p,p1,...,p19;
S e1,...,e19;
V k,l,v,v1,v2,v3,v4,v5;
S s1m,s2m,s3m,s4m,s5m,s6m,s7m,s8m,s9m,s10m,s11m,s12m,s13m,s14m,s15m,s16m,s17m,s18m,s19m;

#do i = 1, `INTLOOPS'
        V p`i'0,...,p`i'19;
        #do j = 1, 19 
                S s`i'`j'm1,...,s`i'`j'm`NUMMASSES';
        #enddo
        S e`i'1,...,e`i'19;
#enddo

#do j = 1, 19 
        S s`j'm1,...,s`j'm`NUMMASSES';
#enddo

I i1,...,i9,j1,...,j9;
I MU,NU,al,be,ro,si,tau,mu,nu,la,ka;

I mu1,...,mu`NUMINDEX', nu1,...,nu`NUMINDEX';
I al1,...,al`NUMINDEX', be1,...,be`NUMINDEX';

#do j = 1, `NUMEXTMOM'

        #ifdef `UPLIMQ`j''
                S pQ`j', pocoQ`j'(:`UPLIMQ`j'');
                #else
                #ifdef `LOWLIMQ`j''
                        S pQ`j', pocoQ`j'(`LOWLIMQ`j'':);
                        #else
                        S pQ`j', pocoQ`j';
                #endif
        #endif
        
#enddo

#do j = 1, `NUMMASSES'

        #ifdef `UPLIMM`j''
                S pM`j', pocoM`j'(:`UPLIMM`j'');
                #else
                #ifdef `LOWLIMM`j''
                        S pM`j', pocoM`j'(`LOWLIMM`j'':);
                        #else
                        S pM`j', pocoM`j';
                #endif
        #endif
        
#enddo

F S,SS,SSS,SSSS, FT1,...,FT`NUMFERMIONTRACES', Slash;
S vmat,vmab,vt,at,vb,ab,ppbt,pmtb,ctt,Htt,ve,ax,sc,ps,nf,nl,nh;
S L,X,exp,[sqrt2],[sqrt(x)],[x],s,g5,g6,g7,ExpZ2;
S x1,...,x9, z2,...,z9,z10;
S xi,x,y,y1,z;
CF ScalProd,ScalProd1,ScalProdM,ScalProdM1,Vec,Den,Den1,Den2;
CF DWt,DWl,VZWW,VZpp,VZWp;
CF Vgh,V3g,Vggs,Dg,Dgh,Dsig;
CF acc,accun,accm;

T FQ,del,dal;

*
* used in MINCER
*

S int1,int2,int3,eq,exp10,exp11,exp20,G311,F321,MSBtoG,y2,y5,j;
S N1,N2,N3,A,B,k1,...,k4;
CF G,poch,po,poinv,ftri;

*
* used in MATAD
*
S dala,diff,test5,test6;
S k5,k6;
S n1,...,n6;
S yy,yy1;
S [p1^2],[M^2+p1^2];
S ggam;
S Mtep;
S S2,OepS2,Oep2S2,T1ep,T1ep2,B4,D3,OepD3,D4,D4ep,D5,D5ep,D6,D6ep,DM,DMep;
S DN,DNep,E3,E3ep;
S intm1,...,intm5,intn1,intt1;
S intbm,intbmbm,intbm1,intbm2,intbn,intbnbn,intbn1,intbn2,intbn3;
S intd4,intd5,intd6,intdm,intdn,inte3,inte4;
S inttbl;
S agam,bgam,cgam,dgam,egam,fgam,ggam,hgam;
CF gm2,gm2norm,gm3,gm3norm;
CF nom,deno,Gam,GGam,iGam,iGGam;
CF BN1;

* AFP
* Exact version
Symbols ep,epp,epp1,epp2,epp3,epp4,epp5,epp6,epp7,epp8,epQ;
CTensor ftensor,dd;
Symbols isum1,isum2,isum3,xpower,n0,n8;

CFunctions rat,acc,den,dena,num,ftriangle;
CFunctions Pochhammer,PochhammerINV,GschemeConstants;
* CFunctions del,ftriangle,hfac;

* two-loop factorized topologies (1-loop)x(1-loop)
S intMxM;
* two-loop integral topologies
S intM00,intMM0,intMMM;
* one-loop topology
S intM0;
* Final topolgy 
S int0;

* eps1m=1/(p1.p1+M^2)^ep
S eps1m,eps2m,eps3m,eps4m,eps5m,eps6m;
* For lower topologies
S epx1,epx2,epx3,epx4,epx5,epx6;

* Zeta values
Symbols z2,z3,z4,z5,z6,z7,z8,z9,zz5,z6z2;

* Master integrals
S miT1,miD6,miD5,miD4,miDN,miE3;
* And its truncation flags
S miT1trunc,miD6trunc,miD5trunc,miD4trunc,miDNtrunc,miE3trunc;
S iGamtrunc,Gamtrunc;

set trunc:miT1trunc,miD6trunc,miD5trunc,miD4trunc,miDNtrunc,miE3trunc,  iGamtrunc,Gamtrunc;

Symbols   [000000], [00000M], [0000M0], [0000MM], 
[000M00], [000M0M], [000MM0], [000MMM], [00M000], 
[00M00M], [00M0M0], [00M0MM], [00MM00], [00MM0M], 
[00MMM0], [00MMMM], [0M0000], [0M000M], [0M00M0], 
[0M00MM], [0M0M00], [0M0M0M], [0M0MM0], [0M0MMM], 
[0MM000], [0MM00M], [0MM0M0], [0MM0MM], [0MMM00], 
[0MMM0M], [0MMMM0], [0MMMMM], [M00000], [M0000M], 
[M000M0], [M000MM], [M00M00], [M00M0M], [M00MM0], 
[M00MMM], [M0M000], [M0M00M], [M0M0M0], [M0M0MM], 
[M0MM00], [M0MM0M], [M0MMM0], [M0MMMM], [MM0000], 
[MM000M], [MM00M0], [MM00MM], [MM0M00], [MM0M0M], 
[MM0MM0], [MM0MMM], [MMM000], [MMM00M], [MMM0M0], 
[MMM0MM], [MMMM00], [MMMM0M], [MMMMM0], [MMMMMM];

CF tad3l,tad2l,tad1l;



PolyRatFun rat;

*
* declarations specific to the problem
*




.global



#procedure ToAuxTopo


* 
*          Our AUX topo is D6
* 
*           
*           _______6_______
*          |\      \      /|
*          |  \       _3/  |
*          |   _\|    /|   |
*         /|1   4 \_/     5|\
*          |      / \      |
*          |    /     \    |
*          |  /         \  |
*          |/______/______\|
*                  2
* 
* 
*         

id tad3l([000000],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

id tad3l([00000M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

id tad3l([0000M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

id tad3l([0000MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4*s5m^n5*s6m^n6;

id tad3l([000M00],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

id tad3l([000M0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3*s4m^n4/p5.p5^n5*s6m^n6;

id tad3l([000MM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3*s4m^n4*s5m^n5/p6.p6^n6;

id tad3l([000MMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3*s4m^n4*s5m^n5*s6m^n6;

id tad3l([00M000],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

id tad3l([00M00M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

id tad3l([00M0M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

id tad3l([00M0MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3/p4.p4^n4*s5m^n5*s6m^n6;

id tad3l([00MM00],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

id tad3l([00MM0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3*s4m^n4/p5.p5^n5*s6m^n6;

id tad3l([00MMM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3*s4m^n4*s5m^n5/p6.p6^n6;

id tad3l([00MMMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3*s4m^n4*s5m^n5*s6m^n6;

id tad3l([0M0000],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

id tad3l([0M000M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

id tad3l([0M00M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

id tad3l([0M00MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3/p4.p4^n4*s5m^n5*s6m^n6;

id tad3l([0M0M00],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

id tad3l([0M0M0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3*s4m^n4/p5.p5^n5*s6m^n6;

id tad3l([0M0MM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3*s4m^n4*s5m^n5/p6.p6^n6;

id tad3l([0M0MMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3*s4m^n4*s5m^n5*s6m^n6;

id tad3l([0MM000],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

id tad3l([0MM00M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

id tad3l([0MM0M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

id tad3l([0MM0MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3/p4.p4^n4*s5m^n5*s6m^n6;

id tad3l([0MMM00],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

id tad3l([0MMM0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3*s4m^n4/p5.p5^n5*s6m^n6;

id tad3l([0MMMM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3*s4m^n4*s5m^n5/p6.p6^n6;

id tad3l([0MMMMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1/p1.p1^n1*s2m^n2*s3m^n3*s4m^n4*s5m^n5*s6m^n6;

id tad3l([M00000],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

id tad3l([M0000M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

id tad3l([M000M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

id tad3l([M000MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3/p4.p4^n4*s5m^n5*s6m^n6;

id tad3l([M00M00],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

id tad3l([M00M0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3*s4m^n4/p5.p5^n5*s6m^n6;

id tad3l([M00MM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3*s4m^n4*s5m^n5/p6.p6^n6;

id tad3l([M00MMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3*s4m^n4*s5m^n5*s6m^n6;

id tad3l([M0M000],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

id tad3l([M0M00M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

id tad3l([M0M0M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

id tad3l([M0M0MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3/p4.p4^n4*s5m^n5*s6m^n6;

id tad3l([M0MM00],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

id tad3l([M0MM0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3*s4m^n4/p5.p5^n5*s6m^n6;

id tad3l([M0MMM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3*s4m^n4*s5m^n5/p6.p6^n6;

id tad3l([M0MMMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1/p2.p2^n2*s3m^n3*s4m^n4*s5m^n5*s6m^n6;

id tad3l([MM0000],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

id tad3l([MM000M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

id tad3l([MM00M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

id tad3l([MM00MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3/p4.p4^n4*s5m^n5*s6m^n6;

id tad3l([MM0M00],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

id tad3l([MM0M0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3*s4m^n4/p5.p5^n5*s6m^n6;

id tad3l([MM0MM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3*s4m^n4*s5m^n5/p6.p6^n6;

id tad3l([MM0MMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2/p3.p3^n3*s4m^n4*s5m^n5*s6m^n6;

id tad3l([MMM000],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3/p4.p4^n4/p5.p5^n5/p6.p6^n6;

id tad3l([MMM00M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3/p4.p4^n4/p5.p5^n5*s6m^n6;

id tad3l([MMM0M0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3/p4.p4^n4*s5m^n5/p6.p6^n6;

id tad3l([MMM0MM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3/p4.p4^n4*s5m^n5*s6m^n6;

id tad3l([MMMM00],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3*s4m^n4/p5.p5^n5/p6.p6^n6;

id tad3l([MMMM0M],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3*s4m^n4/p5.p5^n5*s6m^n6;

id tad3l([MMMMM0],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3*s4m^n4*s5m^n5/p6.p6^n6;

id tad3l([MMMMMM],n1?,n2?,n3?,n4?,n5?,n6?) = 1*s1m^n1*s2m^n2*s3m^n3*s4m^n4*s5m^n5*s6m^n6;

#endprocedure



#procedure Conv2exact()
        
        id po(n?,x?)                = Pochhammer(n,x*ep)*den(x*ep);
        id poinv(n?,x?)             = PochhammerINV(n,x*ep)*num(x*ep);
        
        id nom(x1?,x2?)             = num(x1+x2*ep);
        id deno(x1?,x2?)            = den(x1+x2*ep);

        id nom(x?,y?,z?)            = x + y*num(ep)*num(1+ep*z/y);
*         id nom(x?,y?,z?)            = num(x+y*ep+z*ep^2);    



* id deno(x?,y?) * nom(x?,y?) = 1;
* id deno(0,y?)               = 1/y/ep;
* id nom(x?,y?,z?)            = x + y*nom(0,1)*nom(1,z/y);

* repeat id 1/ep * nom(0,x?) = x;

* id nom(0,y?)         = nom(0,1)*y;
* id deno(x?!{0,1},y?) = deno(1,y/x)/x;
* id nom(x?!{0,1},y?)  = x*nom(1,y/x);
* .sort
        
#endprocedure



#procedure averts(P,in)

        if( count(int`in',1));        
        totensor,nosquare,`P',FQ;

        id FQ(?a) = [sqrt(x)]^nargs_(?a)*FQ(?a);
        id [sqrt(x)]^n?odd_ = 0;
        id,many, [sqrt(x)]*[sqrt(x)] = [x];

        id FQ(?a) = dd_(?a);
        endif;        
        .sort

        if( count(int`in',1));        
        if ( count([x],1) != 0 );
        id [x]^s? =  num(1-ep)*poinv(s+2,-1)*`P'.`P'^(s)/2^s;
        endif;
        endif;        
        .sort
        #call Conv2exact
#endprocedure

#procedure ACCU(TEXT)
*         if ( count(ep,1) != 0);
*         if ( count(acc,1) == 1 ); 
*         id ep^x?*acc(y?) = acc(ep^x*y);
*         else;
*         id ep^x? = acc(ep^x);
*         endif;
*         endif;
*         repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
*         .sort(PolyFun = acc):`TEXT';

        id      nom(x1?,x2?) = num(x1+x2*ep);
        id      deno(x1?,x2?)= den(x1+x2*ep);

        id	num(x?)*den(x?) = 1;
        id	den(x?number_) = 1/x;
        id	num(x?number_) = x;
        id	num(x?) = rat(x,1);
        id	den(x?) = rat(1,x);
        
        .sort:`TEXT';
#endprocedure

* #procedure one10(x,x1,y,y1)

* * 1/(Q.Q + M^2) is "x" and 1/Q.Q is "y";
* * (1/(Q.Q + M^2))^ep is "x1" and (1/Q.Q)^ep is "y1";

*         if ( (count(`x',1) <= 0) && (count(`x1',1) == 0) ) discard;

*         id `x'^k1?*`x1'^k2?*`y'^k3?*`y1'^k4?  = gm2(k1,k2,k3,k4);

* #endprocedure


* #procedure SimpExact

* *         
* *  One-loop tadpole with powers
* *
* *  1/(p.p+m^2)^(1+ep*k2)/p.p^(1+ep*k4)        
* *         
*         id gm2norm(k2?,k4?) =
*         (
*         Gam(1,-1-k4)
*         *Gam(1,1+k2+k4)
*         *iGam(2,-1)
*         *iGam(1,k2)
*         #ifdef 'MINCER'
*                 *ExpZ2
*         #endif
*         );
        
* #endprocedure        

* #procedure simpfin()

* * expand functions w.r.t. ep

* * gm2:

*         id  gm2(k1?,k2?,k3?,k4?) =
*         po(2 - k3    , -k4 - 1)
*         *po(k1 + k3 - 2, k2 + k4 + 1)
*         *poinv(k1     ,k2)
*         *M^(4 - 2*k1 - 2*k3)
* *                  *Exp3(1 + k2 + k4)
*         *gm2norm(k2,k4)
*         ;   

*         id po(1,?a) = 1;
*         id poinv(1,?a) = 1;
*         id po(x1?pos_,0) = fac_(x1-1);
*         id poinv(x1?pos_,0) = 1/(fac_(x1-1));

*         id,many,po(x1?neg0_,x2?) = acc(PO(x1,x2))/x2/ep;
*         id,many,po(x1?,x2?) = acc(PO(x1,x2));

*         repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
*         #call ACCU(simpfin)

*         id,many,poinv(x1?neg0_,x2?) = acc(POINV(x1,x2))*x2*ep;
*         id,many,poinv(x1?,x2?) = acc(POINV(x1,x2));

*         repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
*         #call ACCU(simpfin)

*         id gm2norm(k2?,k4?) =
*         (
*         Gam(1,-1-k4)
*         *Gam(1,1+k2+k4)
*         *iGam(2,-1)
*         *iGam(1,k2)
*         #ifdef 'MINCER'
*                 *ExpZ2
*         #endif
*         );
*         .sort

* ************************************************************

* * gm3:

*         id  gm3(k1?,k2?,k3?,k4?,k5?,k6?) =
*         po(2 - k5    ,      -k6 - 1)
*         *po(k1 + k5 - 2, k2 + k6 + 1)
*         *po(k3 + k5 - 2, k4 + k6 + 1)
*         *po(k1 + k3 + k5 - 4, 2 + k2 +  k4 + k6)
*         *poinv(k1                ,k2)
*         *poinv(k3                ,k4)
*         *poinv(k1 + k3 + 2*k5 - 4,2 + k2 +k4 + 2*k6)
*         *(1 + k2 + k6)
*         *(1 + k4 + k6)
*         *(-1)*nom(1 ,-(2 + k2 + k4 + k6))
*         *(2 + k2 + k4 + k6)
*         /(2 + k2 + k4 + 2*k6)
*         *nom(0,1)*nom(0,1)
*         *M^(8 - 2*k1 - 2*k3 - 2*k5)
*         *gm3norm(k2,k4,k6)
*         ;

*         id po(1,?a) = 1;
*         id poinv(1,?a) = 1;
*         id po(x1?pos_,0) = fac_(x1-1);
*         id poinv(x1?pos_,0) = 1/(fac_(x1-1));

*         id,many,po(x1?neg0_,x2?) = acc(PO(x1,x2))/x2/ep;
*         id,many,po(x1?,x2?) = acc(PO(x1,x2));

*         repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
*         #call ACCU(simpfin)

*         id,many,poinv(x1?neg0_,x2?) = acc(POINV(x1,x2))*x2*ep;
*         id,many,poinv(x1?,x2?) = acc(POINV(x1,x2));

*         repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
*         #call ACCU(simpfin)

*         id gm3norm(x1?,x2?,x3?) = 
*         Gam(-1,2+x1+x2+x3)
*         *Gam(0,1+x1+x3)
*         *Gam(0,1+x2+x3)
*         *Gam(1,-1-x3)
*         *iGam(1,x1 )
*         *iGam(1,x2 )
*         *iGam(0,2+x1+x2+2*x3)
*         *iGam(2,-1)
*         #ifdef 'MINCER'
*                 *ExpZ2^2
*         #endif
*         ;

*         .sort

* ************************************************************

* * l1:

* * Note: The acc's are leftovers when the gamma's are normalized.
* * In that case two gamma's had to be reduced down, so that the
* * first part is one. Actually the first acc can be cancelled
* * against the factors in the definition of the exp's.
* * The last factor comes from the transitition G -> MSbar

*         id G(x1?,y1?,x2?,y2?,n?,s?) = po(x1+x2-s-2,1+y1+y2)
*         *po(2-x1+n-s,-1-y1)*po(2-x2+s,-1-y2)
* *       *acc(1+y1+y2)
*         *nom(1,-2-y1-y2)
*         *poinv(x1,y1)*poinv(x2,y2)*poinv(4-x1-x2+n,-2-y1-y2)
*         *G(1,y1,1,y2,0,0)
*         #ifdef 'MINCER'
*                 *ExpZ2
*         #endif
*         ;
* ***        *Gam(1,1)*Gam(1,-1)^2*iGam(2,-2)
* *
* * the last factor (last line) would be the transitition G -> MSbar
* *


*         id po(1,?a) = 1;
*         id poinv(1,?a) = 1;
*         id po(x1?pos_,0) = fac_(x1-1);
*         id poinv(x1?pos_,0) = 1/(fac_(x1-1));

*         id,many,po(x1?neg0_,x2?) = acc(PO(x1,x2))/x2/ep;
*         id,many,po(x1?,x2?) = acc(PO(x1,x2));

*         repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
*         #call ACCU(simpfin)

*         id,many,poinv(x1?neg0_,x2?) = acc(POINV(x1,x2))*x2*ep;
*         id,many,poinv(x1?,x2?) = acc(POINV(x1,x2));

*         repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
*         #call ACCU(simpfin)


* ***id  G(1,0,1,0,0,0) = 1;
*         id  G(1,0,1,0,0,0) = G(0,0);
*         id  G(1,1,1,1,0,0) = exp11; *1/3
*         id  G(1,0,1,1,0,0) = exp10; *1/2
*         id  G(1,1,1,0,0,0) = exp10; *1/2
*         id  G(1,0,1,2,0,0) = exp20; *1/3
*         id  G(1,2,1,0,0,0) = exp20; *1/3
* ***.sort
*         id  exp11 = G(1,1);
*         id  exp20 = G(2,0);
*         id  exp10 = G(1,0);
*         .sort
*         id G(x1?,x2?)= Gam(1,1+x1+x2)*Gam(1,-1-x1)*Gam(1,-1-x2)
*         *iGam(1,x1)*iGam(1,x2)*iGam(2,-2-x1-x2)
* ***                   *iGam(1,1)*iGam(1,-1)^2*Gam(2,-2)
*         ;


*         #call ACCU(simpfin)


* * Pochhammer functions

*         id po(1,?a) = 1;
*         id poinv(1,?a) = 1;
*         id po(x1?pos_,0) = fac_(x1-1);
*         id poinv(x1?pos_,0) = 1/(fac_(x1-1));

*         id,many,po(x1?neg0_,x2?) = acc(PO(x1,x2))/x2/ep;
*         id,many,po(x1?,x2?) = acc(PO(x1,x2));
*         id,many,poinv(x1?neg0_,x2?) = acc(POINV(x1,x2))*x2*ep;
*         id,many,poinv(x1?,x2?) = acc(POINV(x1,x2));

*         #call ACCU()

* #endprocedure


#procedure tad1l
*
* {tad1l;1;1;0;1; ;(p1:1,1);1;0}
*
        Multiply intM0;
        #call averts(p1,M0)
        #call ACCU{}
        
*         #call one10(s1m|s1m1|1/p1.p1|yy1)
        #call TadpoleM0(s1m,p1,M0,0)        
        .sort

        #call Conv2exact()
        #call DoG
*         #call simpfin()
*         #include expandnomdeno
*         #include redcut
*         #include expepgam
        
#endprocedure        












#procedure partfrac(p1,xxx)
*
* this procedure does the partial fractioning from the term
* 1/p1.p1^a xxx^b where xxx = 1/M^2+p1.p1
*
        id 1/`p1'.`p1' = 1/[p1^2];
        id `xxx' = 1/[M^2+p1^2];  
        ratio [p1^2],[M^2+p1^2],diff;
        id 1/[p1^2] = 1/`p1'.`p1';   
        id 1/[M^2+p1^2] = `xxx';     
        id 1/diff = 1/M^2;
#endprocedure





* #procedure two110(x,x1,y,y1,z,z1)

* * 1/(p1.p1 + M^2) is "x" and 1/(p1.p1 + M^2)^ep  is "x1";
* * 1/(p2.p2 + M^2) is "y" and 1/(p2.p2 + M^2)^ep  is "y1";
* * 1/(p3.p3) is "z" and 1/p3.p3^ep  is "z1";

*         if ( (count(`x',1) <= 0) && (count(`x1',1) == 0) ) discard;
*         if ( (count(`y',1) <= 0) && (count(`y1',1) == 0) ) discard;

*         if ( count(`x1',1)  < count(`y1',1) ) multiply,replace_(`x',`y',`x1',`y1');

*         id `x'^k1?*`x1'^k2?*`y'^k3?*`y1'^k4?*`z'^k5?*`z1'^k6? = gm3(k1,k2,k3,k4,k5,k6);

*         id gm3(k1?,k2?,k3?,k4?,k5?,k6?) =
*         po(2 - k5    ,      -k6 - 1)
*         *po(k1 + k5 - 2, k2 + k6 + 1)
*         *po(k3 + k5 - 2, k4 + k6 + 1)
*         *po(k1 + k3 + k5 - 4, 2 + k2 +  k4 + k6)
*         *poinv(k1                ,k2)
*         *poinv(k3                ,k4)
*         *poinv(k1 + k3 + 2*k5 - 4,2 + k2 +k4 + 2*k6)
*         *(1 + k2 + k6)
*         *(1 + k4 + k6)
*         *(-1)*nom(1 ,-(2 + k2 + k4 + k6))
*         *(2 + k2 + k4 + k6)
*         /(2 + k2 + k4 + 2*k6)
*         *nom(0,1)*nom(0,1)
*         *M^(8 - 2*k1 - 2*k3 - 2*k5)
*         *gm3norm(k2,k4,k6)
*         ;

* #endprocedure





* Exact version
#procedure TadpoleMM0(x,y,z,in,out)

* 1/(p1.p1 + M^2) is "x" and 1/(p1.p1 + M^2)^ep  is "epx";
* 1/(p2.p2 + M^2) is "y" and 1/(p2.p2 + M^2)^ep  is "epy";
* 1/(p3.p3) is "z" and 1/p3.p3^ep  is "epz";
        if(count (int`in',1));
        if ( (count(`x',1) <= 0) && (count(ep`x',1) == 0) ) discard;
        if ( (count(`y',1) <= 0) && (count(ep`y',1) == 0) ) discard;
        
        if ( count(ep`x',1)  < count(ep`y',1) ) multiply,replace_(`x',`y',ep`x',ep`y');
        
        id `x'^k1?*ep`x'^k2?*`y'^k3?*ep`y'^k4?/`z'.`z'^k5?*ep`z'^k6? = gm3(k1,k2,k3,k4,k5,k6);
        
        id gm3(k1?,k2?,k3?,k4?,k5?,k6?) =
        po(2 - k5    ,      -k6 - 1)
        *po(k1 + k5 - 2, k2 + k6 + 1)
        *po(k3 + k5 - 2, k4 + k6 + 1)
        *po(k1 + k3 + k5 - 4, 2 + k2 +  k4 + k6)
        *poinv(k1                ,k2)
        *poinv(k3                ,k4)
        *poinv(k1 + k3 + 2*k5 - 4,2 + k2 +k4 + 2*k6)
        *(1 + k2 + k6)
        *(1 + k4 + k6)
        *(-1)*nom(1 ,-(2 + k2 + k4 + k6))
        *(2 + k2 + k4 + k6)
        /(2 + k2 + k4 + 2*k6)
        *nom(0,1)*nom(0,1)
        *M^(8 - 2*k1 - 2*k3 - 2*k5)
        *gm3norm(k2,k4,k6)
        ;
        Multiply int`out'/int`in';
        endif;
* 
.sort:TadpoleMM0-`in'-1;        
*
        #call Conv2exact
        #call DoG
*         
.sort:TadpoleMM0-`in'-2;        
*         
#endprocedure

#procedure TadpoleM0(x,y,in,out)

* 1/(Q.Q + M^2) is "x" and Q is "y";
* (1/(Q.Q + M^2))^ep is "epx" and (1/Q.Q)^ep is "epy";

        if( count(int`in',1));
        if ( (count(`x',1) <= 0) && (count(ep`x',1) == 0) ) discard;

        id `x'^k1?*ep`x'^k2?/`y'.`y'^k3?*ep`y'^k4?  = gm2(k1,k2,k3,k4);

        id  gm2(k1?,k2?,k3?,k4?) =
        po(2 - k3    , -k4 - 1)*
        po(k1 + k3 - 2, k2 + k4 + 1)*
        poinv(k1     ,k2)*
        M^(4 - 2*k1 - 2*k3)*
*                  Exp3(1 + k2 + k4)*
        gm2norm(k2,k4)
        ;   
        
        Multiply int`out'/int`in';        
        endif;        
#endprocedure

*--#[ DoG :
*
#procedure DoG
*
*	The only objects left are the G(1,x1,1,x2,0,0)
*	which have been written as GschemeConstants(x1,x2)
*
*#$vc = 0;
*Print +f "<1> %t";
id	G(n1?,x1?,n2?,x2?,n3?,n4?) = GschemeConstants(x1,x2)/(1+x1+x2)*
			Pochhammer(n1+n2-n4-2,ep+x1*ep+x2*ep)*
			Pochhammer(1-n1+n3-n4,1-ep-x1*ep)*
			Pochhammer(1-n2+n4,1-ep-x2*ep)*
			PochhammerINV(n1-1,1+x1*ep)*
			PochhammerINV(n2-1,1+x2*ep)*
			PochhammerINV(2-n1-n2+n3,2-2*ep-x1*ep-x2*ep);
*Print +f "<2> %t";
repeat id Pochhammer(n?pos_,x?) = Pochhammer(n-1,x)*num(n-1+x);
repeat id Pochhammer(n?neg_,x?) = Pochhammer(n+1,x)*den(n+x);
repeat id PochhammerINV(n?pos_,x?) = PochhammerINV(n-1,x)*den(n-1+x);
repeat id PochhammerINV(n?neg_,x?) = PochhammerINV(n+1,x)*num(n+x);
id	GschemeConstants(0,n?pos_) = GschemeConstants(n,0);
*	GschemeConstants(0,1)*GschemeConstants(0,2)
*			= GschemeConstants(1,1)*(2*D-6)/(3*D-10)
*			= GschemeConstants(1,1)*(1-2*ep)/(1-3*ep)
id	GschemeConstants(1,1)*GschemeConstants(0,0) =
			GschemeConstants(1,0)*GschemeConstants(2,0)*rat(1-3*ep,1-2*ep);
id	Pochhammer(0,x?) = 1;
id	PochhammerINV(0,x?) = 1;
id	num(x?)*den(x?) = 1;
id	den(x?number_) = 1/x;
id	num(x?number_) = x;
*Print +f "<3> %t";
id	num(x?) = rat(x,1);
*Print +f "<4> %t";
id	den(x?) = rat(1,x);
*Print +f "<5> %t";
*$vc = $vc+1;
*Print +f "<$vc = %$>",$vc;
*
#endprocedure
*
*--#] DoG : 

*--#[ IntOne :
*
#procedure IntOne(p3,p4,Q,in,out)
*
if ( count(int`in',1) );
  if ( ( count(ep`p3',1) == 0 ) && ( count(`p3'.`p3',1) >= 0 ) ) Discard;
  if ( ( count(ep`p4',1) == 0 ) && ( count(`p4'.`p4',1) >= 0 ) ) Discard;
  ToTensor,nosquare,ftensor,`p3';
  if ( count(ftensor,1) == 0 );
	id	int`in'*ep`p3'^x3?*ep`p4'^x4?/`p3'.`p3'^n3?/`p4'.`p4'^n4? =
			int`out'*G(n3,x3,n4,x4,0,0)*`Q'.`Q'^2/`Q'.`Q'^n3/`Q'.`Q'^n4*ep`Q'^x3*ep`Q'^x4*ep`Q';
  elseif ( match(ftensor(i1?)) );
	id	int`in'*ep`p3'^x3?*ep`p4'^x4?/`p3'.`p3'^n3?/`p4'.`p4'^n4?*ftensor(i1?) = int`out'*`Q'(i1)
			*G(n3,x3,n4,x4,1,0)*`Q'.`Q'^2/`Q'.`Q'^n3/`Q'.`Q'^n4*ep`Q'^x3*ep`Q'^x4*ep`Q';
  elseif ( match(ftensor(i1?,i2?)) );
	id	int`in'*ep`p3'^x3?*ep`p4'^x4?/`p3'.`p3'^n3?/`p4'.`p4'^n4?*ftensor(i1?,i2?) =
				int`out'*`Q'.`Q'^2*ep`Q'*ep`Q'^x3*ep`Q'^x4/`Q'.`Q'^n3/`Q'.`Q'^n4*(
			+G(n3,x3,n4,x4,2,0)*`Q'(i1)*`Q'(i2)
			+G(n3,x3,n4,x4,2,1)*d_(i1,i2)*`Q'.`Q'/2);
  else;
	id	int`in'*ep`p3'^x3?*ep`p4'^x4?/`p3'.`p3'^n3?/`p4'.`p4'^n4?*ftensor(?a) = int`out'*ftensor(?a)
			*sum_(isum2,0,integer_(nargs_(?a)/2),G(n3,x3,n4,x4,nargs_(?a),isum2)
				*y^isum2*`Q'.`Q'^isum2/2^isum2)*`Q'.`Q'^2
				*ep`Q'*ep`Q'^x3*ep`Q'^x4/`Q'.`Q'^n3/`Q'.`Q'^n4;
    id  y^isum2?*ftensor(?a) = distrib_(1,2*isum2,del,ftensor,?a);
    tovector,ftensor,`Q';
    id  del(?a) = dd_(?a);
  endif;
  id  P.P = 0;
endif;
*
.sort:IntOne-`in'-1;
*
#call DoG
*
.sort:IntOne-`in'-2;
*
#endprocedure
*
*--#] IntOne : 




#procedure tad2l
*
* {tad2l;3;2;0;1; ;(p1:1,2)(p2:2,1)(p3:2,1);111;110;101;011;100;010;001;000}
*

*
*          /------\
*         /    |   \
*      p1v   p3^    ^p2
*         \    |   /
*          \------/
*             

        #call partfrac{p1|s1m}
        #call partfrac{p2|s2m}
        #call partfrac{p3|s3m}

        if (count(s1m,1,s2m,1,s3m,1) == 0) discard;

* map 010 -> 100
        if ( (count(s1m,1) == 0) && (count(s2m,1) != 0) && (count(s3m,1) == 0) );
        id p1 = -p1;
        id p2 = -p2;
        multiply,replace_(p1,p2,p2,p1,s1m,s2m,s2m,s1m);
        endif;

* map 001 -> 100
        if ( (count(s1m,1) == 0) && (count(s2m,1) == 0) && (count(s3m,1) != 0) );
        id p1 = -p1;
        id p3 = -p3;
        multiply,replace_(p1,p3,p3,p1,s1m,s3m,s3m,s1m);
        endif;

* map 101 -> 110
        if ( (count(s1m,1) != 0) && (count(s2m,1) == 0) && (count(s3m,1) != 0) );
        multiply,replace_(p2,p3,p3,p2,s2m,s3m,s3m,s2m);
        endif;

* map 011 -> 110
        if ( (count(s1m,1) == 0) && (count(s2m,1) != 0) && (count(s3m,1) != 0) );
        id p1 = -p1;
        id p3 = -p3;
        multiply,replace_(p1,p3,p3,p1,s1m,s3m,s3m,s1m);
        endif;

* the rest is taken from topT1 with line 5 replaced by line 3

*
* Now the integration is done:
* (only the momenta p1, p2 and p3 appear)
*

*
* Warning: change direction of line 1 if lines 1 and 2 are massive
*
        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) ); 
        id p1=-p1;
        endif;
        .sort

        #message decompose numerator (M|M|0)

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) ); 
        id p3 = -p1-p2;
        endif;

        #call ACCU{nomgm3 1}

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) ); 
        id p1.p1 = 1/test5 - M^2;
        id p2.p2 = 1/test6 - M^2;
        endif;

        #call ACCU{nomgm3 2}

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) ); 
        id p1.p2 = (p3.p3 - 1/test5 - 1/test6 + 2*M^2)/2;
        endif;

        if (count(s3m,1)>0) id, p3.p3 = 1/s3m - M^2;
        id 1/s3m = p3.p3 + M^2;

        #call ACCU{nomgm3 3}

        #message rec. rel. for the case with 3 massive lines

        #do i=1,1
                #message loop
                if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) ); 
                multiply replace_(test5,s1m,test6,s2m);
                if ( (count(s1m,1)>1) && (count(s2m,1)>=1) && (count(s3m,1)>=1) );
                id s1m^n1?*s2m^n2?*s3m^n5? = 
                s1m^n1*s2m^n2*s3m^n5 * (-1) * 1/3/(n1-1)/M^2 * (
                num(4+3-3*n1 -2*ep)/s1m
                + 2*n2*s2m/s1m*(1/s3m-1/s1m)
                - (n1-1)*(1/s3m-1/s2m)
                );
                redefine i "0";
                endif;
                if ( count(s1m,1) < count(s2m,1) );
                multiply replace_(s1m,s2m,s2m,s1m);
                redefine i "0";
                endif;
                if ( count(s1m,1) < count(s3m,1) );
                multiply replace_(s1m,s3m,s3m,s1m);
                redefine i "0";
                endif;
                if ( count(s2m,1) < count(s3m,1) );
                multiply replace_(s2m,s3m,s3m,s2m);
                redefine i "0";
                endif;
                endif;

                id	num(x?)*den(x?) = 1;
                id	den(x?number_) = 1/x;
                id	num(x?number_) = x;
                id	num(x?) = rat(x,1);
                id	den(x?) = rat(1,x);
                
                .sort                
        #enddo

        #message done
        
        .sort

        #message perform integration

        id 1/s3m = p3.p3 + M^2;

* #include matad.info # time

        #call ACCU{T1}
* #include matad.info # time

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) ); 
        multiply replace_(test5,s1m,test6,s2m);
        Multiply intMM0;
        
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) );
        Multiply intM00;
 
        elseif ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s3m,1)==1) ); 
        id s1m*s2m*s3m = M^2*miT1;
        Multiply int0;
        
        else;
        multiply 1/(1-1);
        endif;
        .sort 
       
* MM0 case
        #call TadpoleMM0(s1m,s2m,p3,MM0,0)
* M00 case        
        #call IntOne(p2,p3,p1,M00,M0)
        #call averts(p1,M0)

* M0 case        
        #call TadpoleM0(s1m,p1,M0,0)

        #call Conv2exact
        #call DoG
        #call subSimple
        #call GammaArgToOne
        
#endprocedure




#procedure bnm2m(TYPE)
*
* Identify M1, M2, ... in the output of BN/BM
*
        if (count(int`TYPE',1) );
*         if (count(intm1,1,intm2,1,intm3,1,intm4,1,intt1,1,intn1,1,
*         intbm,1,intbm1,1,intbm2,1)==0);


        if (count(intm1,1,intm2,1,intm3,1,intm4,1,intt1,1,intn1,1,
        intbm1,1,intbm2,1)==0);
        
        if (match(1/p1.p1/p2.p2*x5*x6)>0); 
        multiply intm1/int`TYPE';
        elseif (match(1/p1.p1/p2.p2/p5.p5*x6)>0); 
        multiply intm2/int`TYPE';
        elseif (match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6)>0); 
        multiply intm3/int`TYPE';
        elseif (match(1/p3.p3/p4.p4*x5*x6)>0); 
        multiply intm4/int`TYPE';
        elseif( (count(p1.p1,1)==0) && (count(p2.p2,1)==0)
        && (count(x3,1)==0) && (count(x4,1)==1)
        && (count(x5,1)==1) && (count(x6,1)==1) );
        multiply intt1/int`TYPE';
        elseif( (count(p1.p1,1)==0) && (count(p2.p2,1)==0)
        && (count(x3,1)==1) && (count(x4,1)==1)
        && (count(x5,1)==1) && (count(x6,1)>0) );
        multiply intn1/int`TYPE';
        endif;
        endif;
        
* BN BM
        endif;
        
#endprocedure        



#procedure symBN1 (p1,p2,x3,x4,p5,x6)
        if( count(intbn1,1));
        if (count(`x3',1) > count(`x4',1))
        multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
        if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) );
        if (count(`x3',1) < count(`x4',1))
        multiply replace_(`x3',`x4',`x4',`x3');
        if (count(`x3',1) < count(`x6',1))
        multiply replace_(`x3',`x6',`x6',`x3');
        if (count(`x4',1) < count(`x6',1))
        multiply replace_(`x4',`x6',`x6',`x4');
        endif;
        endif;
#endprocedure





#procedure redBN1n12to0 (p1,p2,x3,x4,p5,x6)

***#message BN1-n12to0-start

* reduces n1 and n2 to zero, if:
* a. n1<=0, n2>0
* b. n1<=0, n2<0

* a.
*
* reduces n2 to zero if n1<=0 (n3,n4>=1, n5,n6=1)
* result: BN1(n1,n2,n3,n4,n5,n6) with n1<=0, n2=0
*         or BM's.

* sort: n2>=n1

        if( count(intbn1,1) );        
        if ( ( count(`p2'.`p2',1) > count(`p1'.`p1',1) ) &&
        (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)==1) &&
        (count(`p5'.`p5',1)<=-1) )
        multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
        endif;
        
        #call redBN1n6 (`p1',`p2',`x3',`x4',`p5',`x6')
        
        #do i=1,10
                
                if( count(intbn1,1) );        
                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)<=-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)==1) &&
                (count(`p5'.`p5',1)<=-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/( - 2*n5 - n1 + 2*n2 + n4 )
                
                *(
                + `x4'*`x6'^-1 * (  - n4 )
                
                + `p1'.`p1'^-1*`p5'.`p5' * (  - 2*n1 )
                
                + `p1'.`p1'^-1*`x3'^-1 * ( n1 )
                
                + `p1'.`p1'^-1*`x4'^-1 * ( 2*n1 )
                
                + `p1'.`p1'^-1*`x6'^-1 * (  - n1 )
                
                - `p1'.`p1'^-1 * ( 2*n1*M^2 )
                
                + `p2'.`p2'*`x3' * ( 2*n3 )
                
                + `p2'.`p2'*`x4' * ( n4 )
                
                + `p5'.`p5'*`x3' * (  - 2*n3 )
                
                );
                
                redefine i "0";
                
                endif;
                endif;
                #call ACCU(BN1n12to0)
                
        #enddo
        
* #include expandnomdeno

* b.
*
* reduces n2 to zero if n1<=0,n2<0 (n3,n4,n6>=1, n5=1)
* result: BN1(n1,n2,n3,n4,n5,n6) with n1<=0, n2=0
*         or BM's.

        #do i=1,1

* sort: n2>=n1
                if( count(intbn1,1) );
                if ( ( count(`p2'.`p2',1) < count(`p1'.`p1',1) ) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<=-1) )
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');

                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)>0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<=-1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)*deno( + 3/2*4 - n1 - n2 - n3 - n4 - n5 - n6,-3)
                *(
                - `p1'.`p1'^-1*`p2'.`p2'^-1*`x3'^-1 * ( n1*M^2 )

                - `p1'.`p1'^-1*`p2'.`p2'^-1*`x6'^-1 * (  - n1*M^2 )

                - `p2'.`p2'^-1*`p5'.`p5'*`x3' * (  - n3*M^2 )

                + `p2'.`p2'^-1*`x3' * ( 3*n3*M^4 )

                - `p2'.`p2'^-1*`x4'^-1*`x6' * (  - n6*M^2 )

                - `p2'.`p2'^-1*`x4'*`x6'^-1 * (  - n4*M^2 )

                - `p2'.`p2'^-1 *(- 3*n*M^2 + n1*M^2 + 4*n2*M^2 + 3*n3*M^2 + n4*M^2 + 2*
                n5*M^2 + n6*M^2 + 4*M^2 )

                );

                redefine i "0";

                endif;
                endif;
                
                id n=num(4-2*ep);
*                 id acc(x?)*acc(y?)=acc(x*y);
                
*   #include expandnomdeno
                .sort
        #enddo
        
        #call ACCU(BN1n12to0)
        
***#message BN1-n12to0-done
        
#endprocedure

#procedure redBN1n12to1 (p1,p2,x3,x4,p5,x6)

        #do i=1,1
                
* sort: n2>=n1
                if( count(intbn1,1) );
                if ( ( count(`p2'.`p2',1) > count(`p1'.`p1',1) ) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<=-1) )
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                
                if ( (count(`p1'.`p1',1)<=-1) && (count(`p2'.`p2',1)<-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<=-1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)
                *(
                - `x3'^-1 * ( M^-2 )

                - `p1'.`p1'*`p2'.`p2'*`x4' * ( 1/( - 1 + n2)*n4*M^-2 )

                - `p2'.`p2'*`p5'.`p5'*`x4' * (  - 1/( - 1 + n2)*n4*M^-2 )

                + `p2'.`p2'*`x4' * (  - 1/( - 1 + n2)*n4 )

                - `p2'.`p2' * (  - M^-2 + 1/( - 1 + n2)*nom(4,-2)*M^-2 
                - 1/( - 1 + n2)*n4*M^-2 - 2/( - 1 + n2)*n5*M^-2 )

                - `p5'.`p5' * (  - M^-2 )

                );

                redefine i "0";

                endif;
                endif;
                #call ACCU(BN1)
                
        #enddo
        
* #include expandnomdeno
        
***#message BN1-n12to1-done
        
#endprocedure

#procedure redBN1n32 (p1,p2,x3,x4,p5,x6)

***#message BN1-n32-start

* reduce n5 to 1 (if n1=n2=0)
* (take the procedure from redBN1n5.prc)

        #do i=1,1
                
                if( count(intbn1,1) );        
                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2?
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)*M^2*deno(4*n5 - 4 - n1*n5 + n1 + 2*n5 - 2*n5^2,-2*n5+2)
                *(
                + `p1'.`p1'^-1*`p5'.`p5'*`x3'*`x6'^-1 * (  - n1*n3*M^-2 )

                + `p1'.`p1'^-1*`p5'.`p5' * ( n1*n3*M^-2 - n1*n5*M^-2 + n1*M^-2 )

                + `p1'.`p1'^-1*`x4'^-1 * ( n1*n5*M^-2 - n1*M^-2 )

                - `p1'.`p1'^-1 * ( n1*n5 - n1 )

                + `p5'.`p5'*`x3' *(-nom(4,-2)*n3*M^-2+n1*n3*M^-2+2*n3*M^-2+2*n3^2*M^-2 )

                - `p5'.`p5'*`x3'^2 * ( 2*n3 + 2*n3^2 )
                );

                redefine i "0";

                endif;
                endif;
                #call ACCU(BN1n32_0)
                
        #enddo
        
* reduce n3, n4 and n6 to 2, take care if n6=n3-1!!!
* n1=n2=0!
        
* The condition n6=n3-1 is asked for with the following trick; 
* it only works if n5=1!
        
*if ( count(`x3',1) == count(`x6',1,`p5'.`p5',-1) );
*endif;
        
        #do i=1,1
                
                if( count(intbn1,1) );        
                if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) );
                if (count(`x3',1) < count(`x4',1))
                multiply replace_(`x3',`x4',`x4',`x3');
                if (count(`x3',1) < count(`x6',1))
                multiply replace_(`x3',`x6',`x6',`x3');
                if (count(`x4',1) < count(`x6',1))
                multiply replace_(`x4',`x6',`x6',`x4');
                endif;
                
                if ( ( count(`x3',1) > count(`x6',1,`p5'.`p5',-1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(- 6*n3*M^2 + 2*n3^2*M^2 + 4*M^2)
                *(
                - `x3'^-2*`x6' * (  - nom(4,-2)*n6 + n3*n6 + 2*n5*n6 - 2*n6 )
                
                - `x3'^-1*`x4'^-1*`x6' * (  - n3*n6 + 2*n6 )
                
                + `x3'^-1*`x6' * ( n3*n6*M^2 - 2*n6*M^2 )
                
                - `x3'^-1*(4-2*nom(4,-2)*n3+4*nom(4,-2)+2*n3*n5+n3*n6-6*n3+2*n3^2 - 4*n5
                - 2*n6 )
                
                - `p5'.`p5'*`x3'^-1*`x6' * ( n3*n6 - 2*n6 )
                );
*  if (count(`x3',1) < count(`x4',1))
*                   multiply replace_(`x3',`x4',`x4',`x3');
*  if (count(`x3',1) < count(`x6',1))
*                   multiply replace_(`x3',`x6',`x6',`x3');
*  if (count(`x4',1) < count(`x6',1))
*                   multiply replace_(`x4',`x6',`x6',`x4');

                redefine i "0";
                
                endif;
                
                if ( ( count(`x3',1) == count(`x6',1,`p5'.`p5',-1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/( (- 6*n3*M^2 + 2*n3^2*M^2 + 4*M^2) + M^2*n6*(n3-2) )
                *(
                - `x3'^-2*`x6' * (  - nom(4,-2)*n6 + n3*n6 + 2*n5*n6 - 2*n6 )
                
                - `x3'^-1*`x4'^-1*`x6' * (  - n3*n6 + 2*n6 )
                
                - `x3'^-1 * ( 4 - 2*nom(4,-2)*n3 + 4*nom(4,-2) 
                + 2*n3*n5 + n3*n6 - 6*n3 + 2*n3^2 - 4*n5
                - 2*n6 )
                
                - `p5'.`p5'*`x3'^-1*`x6' * ( n3*n6 - 2*n6 )
                );
*  if (count(`x3',1) < count(`x4',1))
*                   multiply replace_(`x3',`x4',`x4',`x3');
*  if (count(`x3',1) < count(`x6',1))
*                   multiply replace_(`x3',`x6',`x6',`x3');
*  if (count(`x4',1) < count(`x6',1))
*                   multiply replace_(`x4',`x6',`x6',`x4');
                
                redefine i "0";

                endif;
                
                if ( ( count(`x3',1) == count(`x6',1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(M^2 * ( 3 - 9/2*n3 + 3/2*n3^2 ) )
                *(
                - 1/`x3'/`x3'*`x6' * (  - nom(4,-2)*n3 + n3^2 )
                
                - 1/`x3'/`x4'*`x6' * ( 2*n3 - n3^2 )
                
                - 1/`x3'*`p5'.`p5'*`x6' * (  - 2*n3 + n3^2 )

                - 1/`x3' * ( 1 - 3/2*nom(4,-2)*n3 + 3*nom(4,-2) - 11/2*n3 + 5/2*n3^2 )
                
                - 1/`x4' * ( 1 - 3/2*n3 + 1/2*n3^2 )
                
                - `p5'.`p5' * (  - 1 + 3/2*n3 - 1/2*n3^2 )
                
                - 1/`x6' * ( 1 + nom(4,-2)*n3 - 2*nom(4,-2) + 5/2*n3 - 3/2*n3^2 )
                
                );
                
*  if (count(`x3',1) < count(`x4',1))
*                   multiply replace_(`x3',`x4',`x4',`x3');
*  if (count(`x3',1) < count(`x6',1))
*                   multiply replace_(`x3',`x6',`x6',`x3');
*  if (count(`x4',1) < count(`x6',1))
*                   multiply replace_(`x4',`x6',`x6',`x4');
                
                redefine i "0";
                
                endif;
                endif;
                
                id nom(4,-2)=num(4-2*ep);
*   id acc(x?)*acc(y?)=acc(x*y);
                #call ACCU(BN1n32_0)
                
        #enddo

***#message BN1-n32-end

#endprocedure

#procedure redBN1n32s (p1,p2,x3,x4,p5,x6)

***#message BN1-n32s-start

* reduce n5 to 1 (if n1=n2=0)
* (take the procedure from redBN1n5.prc)
        
* 13Apr04: 'repeat' -> '#do-#enddo'
***repeat;
        #do ii=1,1
                
                if( count(intbn1,1) );        
                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2?
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)*M^2*deno(4*n5 - 4 - n1*n5 + n1 + 2*n5 - 2*n5^2,-2*n5+2)
                *(
                + `p1'.`p1'^-1*`p5'.`p5'*`x3'*`x6'^-1 * (  - n1*n3*M^-2 )
                
                + `p1'.`p1'^-1*`p5'.`p5' * ( n1*n3*M^-2 - n1*n5*M^-2 + n1*M^-2 )
                
                + `p1'.`p1'^-1*`x4'^-1 * ( n1*n5*M^-2 - n1*M^-2 )
                
                - `p1'.`p1'^-1 * ( n1*n5 - n1 )
                
                + `p5'.`p5'*`x3' *(-nom(4,-2)*n3*M^-2+n1*n3*M^-2+2*n3*M^-2+2*n3^2*M^-2 )
                
                - `p5'.`p5'*`x3'^2 * ( 2*n3 + 2*n3^2 )
                );
                
                redefine j "0";
                redefine ii "0";
                
                endif;
                endif;        
                .sort
                
        #enddo
***endrepeat;
        
        #call ACCU(BN1n32_1)
        
* reduce n3, n4 and n6 to 2, take care if n6=n3-1!!!
* n1=n2=0!
        
* The condition n6=n3-1 is asked for with the following trick; 
* it only works if n5=1!
        
*if ( count(`x3',1) == count(`x6',1,`p5'.`p5',-1) );
*endif;
        
        #do i=1,1

***#message 'i'
                if( count(intbn1,1) );
                if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) );
                if (count(`x3',1) < count(`x4',1))
                multiply replace_(`x3',`x4',`x4',`x3');
                if (count(`x3',1) < count(`x6',1))
                multiply replace_(`x3',`x6',`x6',`x3');
                if (count(`x4',1) < count(`x6',1))
                multiply replace_(`x4',`x6',`x6',`x4');
                endif;
                endif;
                #call ACCU(BN1n32_2)
                
                if( count(intbn1,1) );        
                if ( ( count(`x3',1) > count(`x6',1,`p5'.`p5',-1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(- 6*n3*M^2 + 2*n3^2*M^2 + 4*M^2)
                *(
                - `x3'^-2*`x6' * (  - nom(4,-2)*n6 + n3*n6 + 2*n5*n6 - 2*n6 )
                
                - `x3'^-1*`x4'^-1*`x6' * (  - n3*n6 + 2*n6 )
                
                + `x3'^-1*`x6' * ( n3*n6*M^2 - 2*n6*M^2 )
                
                - `x3'^-1 * ( 4 - 2*nom(4,-2)*n3 + 4*nom(4,-2) + 2*n3*n5 + n3*n6 
                 - 6*n3 + 2*n3^2 - 4*n5
                - 2*n6 )
                
                - `p5'.`p5'*`x3'^-1*`x6' * ( n3*n6 - 2*n6 )
                );
                
                redefine i "0";
                redefine j "0";
                
                endif;
                endif;
                #call ACCU(BN1n32_3)
                
                if( count(intbn1,1) );        
                if ( ( count(`x3',1) == count(`x6',1,`p5'.`p5',-1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/( (- 6*n3*M^2 + 2*n3^2*M^2 + 4*M^2) + M^2*n6*(n3-2) )
                *(
                - `x3'^-2*`x6' * (  - nom(4,-2)*n6 + n3*n6 + 2*n5*n6 - 2*n6 )
                
                - `x3'^-1*`x4'^-1*`x6' * (  - n3*n6 + 2*n6 )
                
                - `x3'^-1 * ( 4 - 2*nom(4,-2)*n3 + 4*nom(4,-2) 
                + 2*n3*n5 + n3*n6 - 6*n3 + 2*n3^2 - 4*n5
                - 2*n6 )
                
                - `p5'.`p5'*`x3'^-1*`x6' * ( n3*n6 - 2*n6 )
                );
                
                redefine i "0";
                redefine j "0";
                
                endif;
                endif;
                
                #call ACCU(BN1n32_4)
                
                if( count(intbn1,1) );        
                if ( ( count(`x3',1) == count(`x6',1) ) &&
                (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)==0) &&
                (count(`x3',1)>2) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)==-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(M^2 * ( 3 - 9/2*n3 + 3/2*n3^2 ) )
                *(
                - 1/`x3'/`x3'*`x6' * (  - nom(4,-2)*n3 + n3^2 )
                
                - 1/`x3'/`x4'*`x6' * ( 2*n3 - n3^2 )
                
                - 1/`x3'*`p5'.`p5'*`x6' * (  - 2*n3 + n3^2 )
                
                - 1/`x3' * ( 1 - 3/2*nom(4,-2)*n3 + 3*nom(4,-2) - 11/2*n3 + 5/2*n3^2 )
                
                - 1/`x4' * ( 1 - 3/2*n3 + 1/2*n3^2 )
                
                - `p5'.`p5' * (  - 1 + 3/2*n3 - 1/2*n3^2 )
                
                - 1/`x6' * ( 1 + nom(4,-2)*n3 - 2*nom(4,-2) + 5/2*n3 - 3/2*n3^2 )
                
                );
                
                redefine i "0";
                redefine j "0";
                
                endif;
                endif;
                
                #call ACCU(BN1n32_5)
                
* added Jul. '98
                
                if( count(intbn1,1) );        
                if ( (count(`p1',1)==0) && (count(`p2',1)==0) );
                if ((count(`x3',1)<=0) || (count(`x4',1)<=0) || (count(`x6',1)<=0)) discard;
                if (count(`x3',1) > count(`x4',1)) multiply replace_(`x3',`x4',`x4',`x3');
                if (count(`x3',1) > count(`x6',1)) multiply replace_(`x3',`x6',`x6',`x3');
                if (count(`x4',1) > count(`x6',1)) multiply replace_(`x4',`x6',`x6',`x4');
                endif;
                endif;
                
                #call ACCU(BN1n32_6)
*   #include expandnomdeno

#enddo

***#message BN1-n32s-done

#endprocedure

#procedure redBN1n34 (p1,p2,x3,x4,p5,x6)

* Use this procedure to reduce n3 and n4 to 1.
* n1=n2=1, n5=n6=1;

* #include expandnomdeno

***#message BN1-n34-start

        #do i=1,1

                if( count(intbn1,1) );        
                if ( ( count(`x3',1) < count(`x4',1) ) &&
                (count(`x3',1)>=1) && (count(`x4',1)>1) && (count(`x6',1)==1) &&
                (count(`p5'.`p5',1)==-1) ) 
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');

                if ( (count(`p1'.`p1',1)==-1) && (count(`p2'.`p2',1)==-1) &&
                (count(`x3',1)>1) && (count(`x4',1)>=1) && (count(`x6',1)==1) &&
                (count(`p5'.`p5',1)==-1) );


                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(n3-1)
                *(
                + `x3'^-1*`x4'^-1 * (  - 3/2*nom(4,-2)*M^-4 
                + n1*M^-4 + n2*M^-4 + n3*M^-4 + n4*M^-4
                + n5*M^-4 + n6*M^-4 - 2*M^-4 )

                - `x3'^-1 * (  - nom(4,-2)*M^-2 + 2*n2*M^-2 
                + n3*M^-2 + n4*M^-2 + n6*M^-2 - 2*M^-2 )

                - `x4'^-1 * ( n3*M^-2 - M^-2 )

                - `p2'.`p2'*`x3'^-1*`x4' * (  - n4*M^-2 )

                + `p2'.`p2'*`x3'^-1 * ( 3/2*nom(4,-2)*M^-4 
                - n1*M^-4 - n2*M^-4 - n3*M^-4 - n4*M^-4 -
                n5*M^-4 - n6*M^-4 + 2*M^-4 )

                - `p5'.`p5' * (  - n3*M^-2 + M^-2 )
                );

                redefine i "0";

                endif;
                endif;
        
                #call ACCU(BN1)

#enddo

* #include expandnomdeno

***#message BN1-n34-done

#endprocedure


#procedure redBN1n5 (p1,p2,x3,x4,p5,x6)
        
        #do i=1,1
                if( count(intbn1,1) );
                if ( (count(`p1'.`p1',1)<=-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(1 - n5)
                *(
                - `x4'^-1 * (  - n5*M^-2 + M^-2 )

                - `p1'.`p1' * ( n5*M^-2 - M^-2 )

                - `p2'.`p2'^-1*`p5'.`p5'*`x4'^-1 * (  - n2*M^-2 )

                - `p2'.`p2'^-1*`p5'.`p5'*`x6'^-1 * ( n2*M^-2 )

                + `p5'.`p5'*`x4' * (  - 2*n4 )

                - `p5'.`p5' * ( nom(4,-2)*M^-2 - n2*M^-2 - 2*n4*M^-2 - n5*M^-2 + M^-2 )
                );

                redefine i "0";

                endif;

                if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)<=-1) &&
                (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
                (count(`p5'.`p5',1)<-1) );

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
                *(-1)/(1 - n5)
                *(
                - `x3'^-1 * (  - n5*M^-2 + M^-2 )

                - `p1'.`p1'^-1*`p5'.`p5'*`x3'^-1 * (  - n1*M^-2 )

                - `p1'.`p1'^-1*`p5'.`p5'*`x6'^-1 * ( n1*M^-2 )

                - `p2'.`p2' * ( n5*M^-2 - M^-2 )

                + `p5'.`p5'*`x3' * (  - 2*n3 )

                - `p5'.`p5' * ( nom(4,-2)*M^-2 - n1*M^-2 - 2*n3*M^-2 - n5*M^-2 + M^-2 )
                );

                redefine i "0";

                endif;
* topBN1        
                endif;
                id nom(4,-2) = num(4-2*ep);

                #call ACCU(BN1)

#enddo

* #include expandnomdeno

#do i=1,1

* #message 'i' n5 (2)

        if( count(intbn1,1) );        
        if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)>=0) &&
        (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
        (count(`p5'.`p5',1)<-1) );

        if (count(`p1'.`p1',1) < count(`p2'.`p2',1))
        multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
        endif;

        if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)==0) &&
        (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
        (count(`p5'.`p5',1)<-1) );

        id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
        * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
        =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
        *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
        *(-1)*M^2*deno(4*n5 - 4 - n1*n5 + n1 + 2*n5 - 2*n5^2,-2*n5+2)
        *(
        + `p1'.`p1'^-1*`p5'.`p5'*`x3'*`x6'^-1 * (  - n1*n3*M^-2 )

        + `p1'.`p1'^-1*`p5'.`p5' * ( n1*n3*M^-2 - n1*n5*M^-2 + n1*M^-2 )

        + `p1'.`p1'^-1*`x4'^-1 * ( n1*n5*M^-2 - n1*M^-2 )

        - `p1'.`p1'^-1 * ( n1*n5 - n1 )

        + `p5'.`p5'*`x3' *(-nom(4,-2)*n3*M^-2 
        + n1*n3*M^-2 + 2*n3*M^-2 + 2*n3^2*M^-2 )

        - `p5'.`p5'*`x3'^2 * ( 2*n3 + 2*n3^2 )
        );

        redefine i "0";

        endif;
        endif;
*   #include expandnomdeno
        .sort
#enddo

***#message BN1-rep-1

if( count(intbn1,1) );
repeat;

        if ( (count(`p1'.`p1',1)<=-1) &&
        (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
        (count(`p5'.`p5',1)<-1) );

        id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
        * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
        =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
        *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
        *(-1)/(1 - n5)
        *(
        - `x4'^-1 * (  - n5*M^-2 + M^-2 )

        - `p1'.`p1' * ( n5*M^-2 - M^-2 )

        - `p2'.`p2'^-1*`p5'.`p5'*`x4'^-1 * (  - n2*M^-2 )

        - `p2'.`p2'^-1*`p5'.`p5'*`x6'^-1 * ( n2*M^-2 )

        + `p5'.`p5'*`x4' * (  - 2*n4 )

        - `p5'.`p5' * ( nom(4,-2)*M^-2 - n2*M^-2 - 2*n4*M^-2 - n5*M^-2 + M^-2 )
        );
        endif;

        if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)<=-1) &&
        (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
        (count(`p5'.`p5',1)<-1) );

        id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
        * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
        =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
        *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
        *(-1)/(1 - n5)
        *(
        - `x3'^-1 * (  - n5*M^-2 + M^-2 )

        - `p1'.`p1'^-1*`p5'.`p5'*`x3'^-1 * (  - n1*M^-2 )

        - `p1'.`p1'^-1*`p5'.`p5'*`x6'^-1 * ( n1*M^-2 )

        - `p2'.`p2' * ( n5*M^-2 - M^-2 )

        + `p5'.`p5'*`x3' * (  - 2*n3 )

        - `p5'.`p5' * ( nom(4,-2)*M^-2 - n1*M^-2 - 2*n3*M^-2 - n5*M^-2 + M^-2 )
        );
        endif;

endrepeat;
endif;
* #include expandnomdeno

***#message BN1-rep-2

if( count(intbn1,1) );
repeat;
        
  if ( (count(`p1'.`p1',1)==0) && (count(`p2'.`p2',1)>=0) &&
     (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
     (count(`p5'.`p5',1)<-1) );

    if (count(`p1'.`p1',1) < count(`p2'.`p2',1))
                   multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
  endif;

  if ( (count(`p1'.`p1',1)>=0) && (count(`p2'.`p2',1)==0) &&
     (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>=1) &&
     (count(`p5'.`p5',1)<-1) );

    id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
      * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
    =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
      *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^n6
      *(-1)*M^2*deno(4*n5 - 4 - n1*n5 + n1 + 2*n5 - 2*n5^2,-2*n5+2)
      *(
       + `p1'.`p1'^-1*`p5'.`p5'*`x3'*`x6'^-1 * (  - n1*n3*M^-2 )

       + `p1'.`p1'^-1*`p5'.`p5' * ( n1*n3*M^-2 - n1*n5*M^-2 + n1*M^-2 )

       + `p1'.`p1'^-1*`x4'^-1 * ( n1*n5*M^-2 - n1*M^-2 )

       - `p1'.`p1'^-1 * ( n1*n5 - n1 )

       + `p5'.`p5'*`x3' *(-nom(4,-2)*n3*M^-2 
                  + n1*n3*M^-2 + 2*n3*M^-2 + 2*n3^2*M^-2 )

       - `p5'.`p5'*`x3'^2 * ( 2*n3 + 2*n3^2 )
       );
  endif;

endrepeat;
endif;

***#message BN1-rep-done

#endprocedure

#procedure redBN1n6 (p1,p2,x3,x4,p5,x6)
        
        if( count(intbn1,1) );        
        repeat;
                if ( (count(`x3',1)>=1) && (count(`x4',1)>=1) && (count(`x6',1)>1) &&
                (count(`p5'.`p5',1)<=-1) );
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * 1/`p5'.`p5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * 1/`p5'.`p5'^n5  * `x6'^(n6-1)
                *1/(n6-1)
                *(
                - n3*`x3' - n4*`x4'
                - 1/M^2*nom(6-n1-n2-n3-n4-n5-(n6-1),-3)
                );
                endif;
        endrepeat;
        endif;

#endprocedure





#procedure reduceBN1
        
        #call redBN1n5(p1,p2,x3,x4,p5,x6)
        .sort
***#call symBN1(p1,p2,x3,x4,p5,x6)
        
        #call redBN1n12to1(p1,p2,x3,x4,p5,x6)
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
        #call redBN1n6 (p1,p2,x3,x4,p5,x6)
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
        #call redBN1n34 (p1,p2,x3,x4,p5,x6)
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
* #include expandnomdeno
        
        #call redBN1n12to0 (p1,p2,x3,x4,p5,x6)
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
* set massless tadpoles to zero:
        
        if( count(intbn1,1) );        
        if ( (count(x3,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x4,1)<=0) &&
        ( (count(p1.p1,1)>=0) || (count(p2.p2,1)>=0) || (count(p5.p5,1)>=0) )
        ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) && (count(p5.p5,1)<0) &&
        ( (count(x3,1)<=0)    || (count(x4,1)<=0) ||
        (count(x6,1)<=0) )
        ) discard;
        endif;        
        .sort
        
* #include expandnomdeno
        
        #do j=1,1
                #call redBN1n32s (p1,p2,x3,x4,p5,x6)
        #enddo
        
        #call redBN1n32 (p1,p2,x3,x4,p5,x6)
        .sort
        
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
* now treat the integrals BN1(0,0,2,2,1,1), BN1(0,0,2,2,1,2)
* and BN1(0,0,2,1,1,1) separate:
        
* BN1(0,0,2,2,1,1) = -1/3/M^2*BN1(0,0,2,1,0,2);
        
        if( count(intbn1,1) );
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) && (count(x3,1)==2) && 
        (count(x4,1)==2) && (count(p5.p5,1)==-1) && (count(x6,1)==1) )
        id x3^2*x4^2/p5.p5*x6 = -1/3/M^2 * x3^2*x4*x6^2;
        endif;
        .sort
        
*   BN(0,0,2,2,1,2) =
*       + 7/6*BN1(0,0,2,1,0,2)*n*M^-4 - 38/9*BN1(0,0,2,1,0,2)*M^-4 
*       + 19/3*BN1(0,0,2,1,1,1)*n*M^-4 - BN1(0,0,2,1,1,1)*n^2*M^-4 
*       - 10*BN1(0,0,2,1,1,1)*M^-4 + 4/3*BN1(0,0,3,1,0,2)*M^-2;
        
        if( count(intbn1,1) );
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) && (count(x3,1)==2) && 
        (count(x4,1)==2) && (count(p5.p5,1)==-1) && (count(x6,1)==2) )
        id x3^2*x4^2/p5.p5*x6^2 =
        + 7/6* x3^2*x4*x6^2 *n*M^-4      - 38/9* x3^2*x4*x6^2 *M^-4 
        + 19/3* x3^2*x4/p5.p5*x6 *n*M^-4 - x3^2*x4/p5.p5*x6 *n^2*M^-4 
        - 10* x3^2*x4/p5.p5*x6 *M^-4     + 4/3* x3^3*x4*x6^2 *M^-2
        ;
        endif;
        .sort
        
*   BN(0,0,2,1,1,1) =
*      -1/3/M^2*(3/2*n-4)*BN1(0,0,1,1,1,1);
        
        if( count(intbn1,1) );
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) && (count(x3,1)==2) && 
        (count(x4,1)==1) && (count(p5.p5,1)==-1) && (count(x6,1)==1) )
        id x3^2*x4/p5.p5*x6 =
        -1/3/M^2*(3/2*nom(4,-2)-4)* x3*x4/p5.p5*x6
        ;
        endif;
        
        id nom(4,-2)=num(4-2*ep);
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
#endprocedure        





#procedure topbn1
*
* this is topbn1
*
        #-
        #message this is topbn1
        
        #message numerator
        
        if( count(intbn1,1) );        
        id  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
        id  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
        id  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
        endif;        
        #call ACCU(BN1 1)
        
        if( count(intbn1,1) );                
        id  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );
        id  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
        endif;        
        #call ACCU(BN1 2)
        
        if( count(intbn1,1) );                
        id  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
        id  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
        endif;        
        #call ACCU(BN1 3)
        
        if( count(intbn1,1) );                
        id  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
        id  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
        endif;        
        #call ACCU(BN1 4)
        
        if( count(intbn1,1) );        
        id  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
        id  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
        endif;        
        #call ACCU(BN1 5)
        
        if( count(intbn1,1) );        
        id  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
        id  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN1 6)
        
        if( count(intbn1,1) );        
        id  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
        id  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN1 7)
        
*
* Warning!
*
        if( count(intbn1,1) );
        id  1/x5 = M^2 + p5.p5;
        
        id  p3.p3 = 1/x3 - M^2;
        endif;
        
        #call ACCU(BN1 6)
        if( count(intbn1,1) );        
        id  p4.p4 = 1/x4 - M^2;
        id  p6.p6 = 1/x6 - M^2;
        endif;
        
        #call ACCU(BN1 7)
        
        #call symBN1(p1,p2,x3,x4,p5,x6)
        .sort
        
        if( count(intbn1,1) );        
        id x3^n3?neg_=(p3.p3+M^2)^-n3;
        id x4^n4?neg_=(p4.p4+M^2)^-n4;
        id x6^n6?neg_=(p6.p6+M^2)^-n6;
        
        
        if ( (count(x3,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x4,1)<=0) && 
        ( (count(p1.p1,1)>=0) || (count(p2.p2,1)>=0) || (count(p5.p5,1)>=0) )
        ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) && (count(p5.p5,1)<0) &&
        ( (count(x3,1)<=0)    || (count(x4,1)<=0) || 
        (count(x6,1)<=0) )
        ) discard;
        endif;        
        .sort
        
        #message do recursion
        
        #call reduceBN1
        
        if( count(intbn1,1) );        
        if ( (count(x3,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(p5.p5,1)>=0) ) discard;
        if ( (count(x3,1)<=0) && (count(x4,1)<=0) && 
        ( (count(p1.p1,1)>=0) || (count(p2.p2,1)>=0) || (count(p5.p5,1)>=0) )
        ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) && (count(p5.p5,1)<0) &&
        ( (count(x3,1)<=0)    || (count(x4,1)<=0) || 
        (count(x6,1)<=0) )
        ) discard;
        
        if (count(x6,1)==0);
        id p6=p1+p3;
        id p1=-p1;
        multiply replace_(p5,p4,p4,p6,p3,p5,x3,x5,x4,x6);
        multiply intbm1/intbn1;
        elseif (count(x3,1)==0);
        id p3=p5-p2;
        id p1=-p1;
        multiply replace_(p5,p4,p4,p6,p6,p5,p2,p3,x4,x6,x6,x5);
        multiply intbm1/intbn1;
        elseif (count(x4,1)==0);
        id p4=p6-p2;
        multiply replace_(p5,p4,p3,p5,p1,p3,x3,x5);
        multiply intbm1/intbn1;
        elseif (count(p5.p5,1)>=0);
        id p5=p4-p1;
        multiply replace_(p3,p5,p1,p3,p2,p1,x3,x5);
        multiply intbm/intbn1;
        elseif ( (count(x3,1)!=0) && (count(x4,1)!=0) && (count(x6,1)!=0) );
        id 1/p1.p1^n1? /p2.p2^n2? * x3^n3? * x4^n4? /p5.p5^n5? * x6^n6? =
        BN1(n1,n2,n3,n4,n5,n6);        
        endif;
* topBN1        
        endif;        
.sort

* insert the expansion for BN1(0,0,1,1,1,1) and BN1(1,1,1,1,1,1).

* id BN1(0,0,1,1,1,1) = +  M^4*(
*       + ep^-3 + ep^-2 * ( 15/4 ) + ep^-1 * ( 65/8 + 12/8*z2 )
*       + 81/4*S2 - z3 + 135/16 + 90/16 *z2
*       + ep*OepS2
*       + ep^2*Oep2S2
*                            )
* #ifdef 'MINCER'
*       *(1 + z2*ep^2/2 + 1/8*z2^2*ep^4 + 1/48*z2^3*ep^6)^3;
* #endif
* ;

* id BN1(1,1,1,1,1,1) = + (2*z3/ep + D3) + ep*OepD3;
* .sort

#call ACCU();

#message - done
        
#endprocedure        





#procedure dorec3l
*
* dorec3l -> master procedure for three loop tadpole recurrence relations
* 
************************************************************ 

* Note:
* - Treat first d6, d5, ...; then bn, ... and at last bm, ...
* - MINCER calculates the diagrams originally in the G-scheme
*   and changes at the end to the MSbar-scheme. In addition 
*   a multiplication of exp(z2*\ep^2/2) per loop is done.
* - The topologies BN1 and BN2 are reduced to topologies BM, BM1 
*   and simpler functions BN1(...) (of course topology BN1 only, 
*   this should be changed in future), so terms proportional to
*   intbm, intbm1 and BN1(...) are present inside intbn1 and intbn2
*   which we have to be set equal to zero, in order not to count anything
*   twice. (BN1(...) are added to 'diarest'; see below.) 

************************************************************ 

*         multiply replace_(intbn,intbnbn,intbm,intbmbm);
*         .sort

*         Print+s;
* .end        
        
        #do type = {d6|d5|d4|dm|dn|e4|e3}

*   G dia`type' = dia;

*   if (count(int`type',1) == 0) discard;
*   id int`type' = 1;
*   .sort

***  #if (termsin(dia`type')!=0)

                #message Recursion of type `type'

*     #include matad.info # time 
*     #include matad.info # print

                #call top`type'
                Print+s;        
                .sort        
*     #include matad.info # time
*     #include matad.info # print

*     G diatmp = dia + dia`type';
                
*     multiply replace_(intbn,intbnbn,intbm,intbmbm);
*     id int`type' = 0;
*     .sort

*   drop dia`type';
*   .sort 
*   delete storage;
*   .sort
*   .store

*   g dia = diatmp;
*   .sort
*   delete storage;
*   .sort
*   .store

#enddo

* #do type = {bnbn|bn1|bn2|bn3|bmbm|bm1|bm2|m1|m2|m3|m4|m5|t1|n1}


#do type = {bn|bn1|bn2|bn3|bm|bm1|bm2|m1|m2|m3|m4|m5|t1|n1}

        
*   G dia`type' = dia;

*   if (count(int`type',1)==0) discard;
*   id int`type' = 1;

*   .sort

* ***  #if (termsin(dia`type')!=0)

    #message Recursion of type `type'
    multiply replace_(s1m,x1,s2m,x2,s3m,x3,s4m,x4,s5m,x5,s6m,x6,s7m,x7,s8m,x8);

* *     #include matad.info # time 
* *     #include matad.info # print

    #call top`type'
    #call bnm2m(`type')

* *     #include matad.info # time
* *     #include matad.info # print

*     G diatmp = dia + dia`type';

*         multiply replace_(intbn,intbnbn,intbm,intbmbm);
*     id int`type' = 0;
*     .sort

*   drop dia`type';
*   .sort
*   delete storage;
*   .sort
*   .store
*   g dia = diatmp;
*   .sort
*   delete storage;
*   .sort
*   .store

#enddo

* ************************************************************

*   g diatmp = dia;
*   .sort
*   delete storage;
*   .sort
*   .store
*   g dia = diatmp;
*   .sort
*   delete storage;
*   .sort

* #message final manipulations

* * #include matad.info # time
* * #include matad.info # print

* #message simplify

* #call simpfin
* #include redcut

* #call ACCU()

* * use "expandnomdeno" because of treatment of 1/ep poles

* #include expandnomdeno2

* * expand ep

* #include expepgam

* * #include matad.info # time
* * #include matad.info # print

* id acc(x?)=x;
* .sort 
#endprocedure






************************************************************

#procedure mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)
        
* discard massless tadpoles
        if ( count(intd6,1) );
        if ( (count(`s1m',1)<=0)&&(count(`s2m',1)<=0)&&(count(`s4m',1)<=0) ) discard;
        if ( (count(`s1m',1)<=0)&&(count(`s3m',1)<=0)&&(count(`s6m',1)<=0) ) discard;
        if ( (count(`s2m',1)<=0)&&(count(`s3m',1)<=0)&&(count(`s5m',1)<=0) ) discard;
        if ( (count(`s4m',1)<=0)&&(count(`s5m',1)<=0)&&(count(`s6m',1)<=0) ) discard;
        endif;        
        .sort

#endprocedure

#procedure redD6n5(s1m,s2m,s3m,s4m,s5m,s6m)

* reduce n5 from >1 to =1

        if (match(`s1m'*`s2m'*`s3m'*`s4m'*`s5m'^2*`s6m') > 0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?*`s6m'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5*`s6m'^n6 * 
        (-1)/4/(n5-1)/M^2/`s5m' * (
        n-4*(n5-1)
        -3*( n3*`s3m'*(1/`s5m'-1/`s6m') + n2*`s2m'*(1/`s5m'-1/`s4m') )
        +n2*`s2m'*(1/`s3m'-1/`s1m') 
        +n3*`s3m'*(1/`s2m'-1/`s1m') 
        +(n5-1)*`s5m'*(1/`s3m'-1/`s6m'+1/`s2m'-1/`s4m')
        );

* sort: n5 >= n1,n2,n3,n4,n6

        if (count(`s1m',1) > count(`s5m',1)) 
        multiply replace_(`s1m',`s3m',`s3m',`s1m',`s4m',`s5m',`s5m',`s4m');
        if (count(`s5m',1) < count(`s6m',1)) 
        multiply replace_(`s5m',`s6m',`s6m',`s5m',`s1m',`s2m',`s2m',`s1m');
        if (count(`s5m',1) < count(`s4m',1)) 
        multiply replace_(`s5m',`s4m',`s4m',`s5m',`s1m',`s3m',`s3m',`s1m');
        if (count(`s5m',1) < count(`s3m',1)) 
        multiply replace_(`s5m',`s3m',`s3m',`s5m',`s1m',`s4m',`s4m',`s1m');
        if (count(`s5m',1) < count(`s2m',1)) 
        multiply replace_(`s5m',`s2m',`s2m',`s5m',`s1m',`s6m',`s6m',`s1m');

        redefine i "0";

        endif;

#endprocedure

************************************************************

* treat the scalar products

#procedure topd6
*
* this is topd6
*

        #message this is topd6

        #message numerator

* AFP        
        if ( count(intd6,1) );
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p5.p5 = 1/s5m - M^2;
        id p6.p6 = 1/s6m - M^2;
        endif;
        #call ACCU(D6 0)
        #call mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)

        if ( count(intd6,1) );        
        id p2 = p1+p3;
        id p4 = p1+p6;
        id p5 = p6-p3;
        endif;        
        #call ACCU(D6 1)

        if ( count(intd6,1) );        
        id p1.p3 = 1/2 * (   1/s2m - 1/s1m - 1/s3m + M^2 );
        id p1.p6 = 1/2 * (   1/s4m - 1/s1m - 1/s6m + M^2 );
        id p3.p6 = 1/2 * ( - 1/s5m + 1/s3m + 1/s6m - M^2 );
        endif;        
        #call ACCU(D6 2)
        
        #call mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)

        if ( count(intd6,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p5.p5 = 1/s5m - M^2;
        id p6.p6 = 1/s6m - M^2;
        endif;        
        #call ACCU(D6 3)
        
        #call mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)
        
        #call ACCU(D6)

* #include matad.info # time
* #include matad.info # print
        
************************************************************
        
* do recursion
        
        #message do recursion
        
* sort: n5 >= n1,n2,n3,n4,n6
        if ( count(intd6,1) );
        if (count(s1m,1) > count(s5m,1)) 
        multiply replace_(s1m,s3m,s3m,s1m,s4m,s5m,s5m,s4m);
        if (count(s5m,1) < count(s6m,1)) 
        multiply replace_(s5m,s6m,s6m,s5m,s1m,s2m,s2m,s1m);
        if (count(s5m,1) < count(s4m,1)) 
        multiply replace_(s5m,s4m,s4m,s5m,s1m,s3m,s3m,s1m);
        if (count(s5m,1) < count(s3m,1)) 
        multiply replace_(s5m,s3m,s3m,s5m,s1m,s4m,s4m,s1m);
        if (count(s5m,1) < count(s2m,1)) 
        multiply replace_(s5m,s2m,s2m,s5m,s1m,s6m,s6m,s1m);
        endif;        

        #do i=1,1
                #call redD6n5(s1m,s2m,s3m,s4m,s5m,s6m)
                .sort
        #enddo
        
        #call mltadD6(s1m,s2m,s3m,s4m,s5m,s6m)
        
        #call ACCU(D6)
        
************************************************************
        
* identify simpler integrals
        if ( count(intd6,1) );
        if ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2; 
        id p1=p4-p6;
        id p2=-p2;
        id p4=-p4;
        multiply replace_(p2,p1,p3,p4,p4,p3,p5,p2,p6,p5,
	s2m,s1m,s3m,s4m,s4m,s3m,s5m,s2m,s6m,s5m);
        multiply, intd5/intd6;
        elseif ( (count(s2m,1)<=0) );
        id 1/s2m=p2.p2+M^2; 
        id p2=p1+p3;
        id p1=-p1;
        id p3=-p3;
        multiply replace_(p1,p3,p3,p5,p4,p1,p5,p4,p6,p2,
	s1m,s3m,s3m,s5m,s4m,s1m,s5m,s4m,s6m,s2m);
        multiply, intd5/intd6;
        elseif ( (count(s3m,1)<=0) );
        id 1/s3m=p3.p3+M^2; 
        id p3=p6-p5;
        id p1=-p1;
        id p2=-p2;
        id p4=-p4;
        id p6=-p6;
        multiply replace_(p2,p4,p4,p2,p6,p3,
        s2m,s4m,s4m,s2m,s6m,s3m);
        multiply, intd5/intd6;
        elseif ( (count(s4m,1)<=0) );
        id 1/s4m=p4.p4+M^2; 
        id p4=p2+p5;
        id p2=-p2;
        id p3=-p3;
        multiply replace_(p1,p3,p2,p1,p3,p2,p5,p4,p6,p5,
        s1m,s3m,s2m,s1m,s3m,s2m,s5m,s4m,s6m,s5m);
        multiply, intd5/intd6;
        elseif ( (count(s5m,1)<=0) );
        id 1/s5m=p5.p5+M^2; 
        id p5=p6-p3;
        id p1=-p1;
        id p2=-p2;
        multiply replace_(p1,p2,p2,p3,p3,p1,p4,p5,p6,p4,
        s1m,s2m,s2m,s3m,s3m,s1m,s4m,s5m,s6m,s4m);
        multiply, intd5/intd6;
        elseif ( (count(s6m,1)<=0) );
        id 1/s6m=p6.p6+M^2; 
        id p6=p4-p1;
        multiply, intd5/intd6;
        elseif ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s3m,1)==1) && 
        (count(s4m,1)==1) && (count(s5m,1)==1) && (count(s6m,1)==1) );
        id s1m*s2m*s3m*s4m*s5m*s6m = miD6;
        Multiply int0/intd6;
        else;
        multiply 1/(1-1);
        endif;
        endif;
        
* #include expandnomdeno

        #message - done
        Print+s;
        .sort        
#endprocedure        



************************************************************

#procedure mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

* discard massless tadpoles
        if ( count(intd5,1) );
        if ( (count(`s1m',1)<=0) && (count(`s3m',1)<=0) ) discard;
        if ( (count(`s4m',1)<=0) && (count(`s5m',1)<=0) ) discard;
        endif;        
        .sort
        
#endprocedure

#procedure redD5n6p(s1m,s2m,s3m,s4m,s5m,p6)

* reduce n6 from >0 to =0
        if ( count(intd5,1) );        
        if (match(`s1m'*`s2m'*`s3m'*`s4m'*`s5m'/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 
        
        (-1)*deno(4-2*n6-n3-n1,-2)*(
        + n1*`s1m'*`s4m'^-1
        + n3*`s3m'*`s5m'^-1
        - `p6'.`p6'*n1*`s1m'
        - `p6'.`p6'*n3*`s3m'
        );
        
        endif;
* topd5                
        endif;        
#endprocedure

#procedure redD5n6m(s1m,s2m,s3m,s4m,s5m,p6)

* reduce n6 from <0 to =0
        if ( count(intd5,1) );        
        if (match(`s1m'*`s2m'*`s3m'*`s4m'*`s5m'*`p6'.`p6')>0);
        
        if (count(`s3m',1)==1);
        
        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 
        
        (-1)*deno(4-2*n3-n6-n1,-2)*(
        + n1*M^2*`s1m'                 
        + n1*`s1m'*`s2m'^-1            
        - n1*`s1m'*`s3m'^-1            
        + 2*n3*M^2*`s3m'               
        - `p6'.`p6'^-1*n6*`s3m'^-1     
        + `p6'.`p6'^-1*n6*`s5m'^-1     
        );
        
* sort: n3 >= n4,n5,n1
        
        repeat;
                if (count(`s3m',1) < count(`s4m',1)) 
                multiply replace_(`s3m',`s4m',`s4m',`s3m',`s1m',`s5m',`s5m',`s1m');
                if (count(`s3m',1) < count(`s5m',1)) 
                multiply replace_(`s3m',`s5m',`s5m',`s3m',`s1m',`s4m',`s4m',`s1m');
                if (count(`s3m',1) < count(`s1m',1)) 
                multiply replace_(`s3m',`s1m',`s1m',`s3m',`s5m',`s4m',`s4m',`s5m');
        endrepeat;
        
        redefine i "0";
        
        endif;

        endif;


        if (match(`s1m'*`s2m'*`s3m'*`s4m'*`s5m'*`p6'.`p6')>0);
        
* 29Oct04: added
        if (count(`s3m',1)!=1);
        
        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 
        
        (-1)/(n3-1)*(
        + `p6'.`p6'^-1*n2*M^2*`s2m'*`s3m'^-1     
        + `p6'.`p6'^-1*n2*`s2m'*`s3m'^-1*`s4m'^-1 
        - `p6'.`p6'^-1*n2*`s2m'*`s3m'^-1*`s5m'^-1
        - `p6'.`p6'^-1*n2*`s3m'^-1
        + 2*`p6'.`p6'^-1*n5*M^2*`s3m'^-1*`s5m'
        - 2*`p6'.`p6'^-1*n5*`s3m'^-1
        + 2*`p6'.`p6'^-1*M^2*(n3-1)
        - `p6'.`p6'^-1*(n3-1)*`s3m'^-1
        - `p6'.`p6'^-1*(n3-1)*`s5m'^-1
        + `p6'.`p6'^-1*`s3m'^-1*n
        );
        
* 29Oct04: added
        endif;
        
* 29Oct04: added
        redefine i "0";
        
        endif;
* topd5        
        endif;        
#endprocedure


#procedure redD5n2(s1m,s2m,s3m,s4m,s5m,p6)

* reduce n2 from >1 to =1
* applied for n6=0

        if ( count(intd5,1) );        
        if ( (match(`s1m'*`s2m'^2*`s3m'*`s4m'*`s5m')>0) && (count(`p6',1)==0) );
        
        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 
        
        (-1)/(n2-1)*(
        - n1*M^-2*`s2m'^-1
        + n1*`s1m'*`s2m'^-1
        - n3*M^-2*`s2m'^-1
        + n3*`s2m'^-1*`s3m'
        - n4*M^-2*`s2m'^-1
        + n4*`s2m'^-1*`s4m'
        - n5*M^-2*`s2m'^-1
        + n5*`s2m'^-1*`s5m'
        - n6*M^-2*`s2m'^-1
        - M^-2*(n2-1)*`s2m'^-1
        + 3/2*M^-2*`s2m'^-1*n
        );
        
        redefine i "0";
        
        endif;
* topd5        
        endif;
#endprocedure


#procedure redD5n1345(s1m,s2m,s3m,s4m,s5m,p6)

* reduce n1,n3,n4 and n5 from >1 to =1
* requirements: n2=1 and n6=0 

        if ( count(intd5,1) );        
        if ( (match(`s1m'^2*`s2m'*`s3m'*`s4m'*`s5m')>0) && 
        (count(`p6',1)==0) && (count(`s2m',1)==1) );

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?*`s4m'^n4?*`s5m'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3*`s4m'^n4*`s5m'^n5/`p6'.`p6'^n6 * 

        (-1)/(n1-1)*(
        - 2/3*n3*M^-2*`s1m'^-2*`s3m'
        + 2/3*n3*M^-2*`s1m'^-1*`s2m'^-1*`s3m'
        - 1/3*n6*M^-2*`s1m'^-1
        - M^-2*(n1-1)*`s1m'^-1
        - 1/3*M^-2*(n1-1)*`s2m'^-1
        + 1/3*M^-2*(n1-1)*`s3m'^-1
        + 1/3*M^-2*`s1m'^-1*n
        - 2/3*`p6'.`p6'^-1*n6*M^-2*`s1m'^-2
        + 1/3*`p6'.`p6'^-1*n6*M^-2*`s1m'^-1*`s3m'^-1
        + 2/3*`p6'.`p6'^-1*n6*M^-2*`s1m'^-1*`s4m'^-1
        - 1/3*`p6'.`p6'^-1*n6*M^-2*`s1m'^-1*`s5m'^-1
        );

* sort: n1 >= n3,n4,n5

        repeat;
                if (count(`s1m',1) < count(`s3m',1)) 
                multiply replace_(`s1m',`s3m',`s3m',`s1m',`s4m',`s5m',`s5m',`s4m');
                if (count(`s1m',1) < count(`s4m',1)) 
                multiply replace_(`s1m',`s4m',`s4m',`s1m',`s3m',`s5m',`s5m',`s3m');
                if (count(`s1m',1) < count(`s5m',1)) 
                multiply replace_(`s1m',`s5m',`s5m',`s1m',`s3m',`s4m',`s4m',`s3m');
        endrepeat;

        if ( (count(`s2m',1)==0) && (count(`p6',1)==0) );
        if (count(`s1m',1)<count(`s3m',1)) multiply replace_(`s1m',`s3m',`s3m',`s1m');
        if (count(`s1m',1)<count(`s4m',1)) multiply replace_(`s1m',`s4m',`s4m',`s1m');
        if (count(`s1m',1)<count(`s5m',1)) multiply replace_(`s1m',`s5m',`s5m',`s1m');
        if (count(`s3m',1)<count(`s4m',1)) multiply replace_(`s3m',`s4m',`s4m',`s3m');
        if (count(`s3m',1)<count(`s5m',1)) multiply replace_(`s3m',`s5m',`s5m',`s3m');
        if (count(`s4m',1)<count(`s5m',1)) multiply replace_(`s4m',`s5m',`s5m',`s4m');
        endif;

        if ( (count(`s3m',1)==0) && (count(`p6',1)==0) );
        if (count(`s2m',1)<count(`s4m',1)) multiply replace_(`s2m',`s4m',`s4m',`s2m');
        if (count(`s2m',1)<count(`s5m',1)) multiply replace_(`s2m',`s5m',`s5m',`s2m');
        if (count(`s4m',1)<count(`s5m',1)) multiply replace_(`s4m',`s5m',`s5m',`s4m');
        endif;

        redefine i "0";

        endif;
* topd5  
        endif;
#endprocedure

#procedure topd5
*
* this is topd5
*
        #message this is topd5

************************************************************

* treat the scalar products

        #message numerator

        if ( count(intd5,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p5.p5 = 1/s5m - M^2;
        endif;        
        #call ACCU(D5 0)

        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)
        
        if ( count(intd5,1) );        
        id p2 = p1+p3;
        id p4 = p1+p6;
        id p5 = p6-p3;
        endif;        
        #call ACCU(D5 1)

        if ( count(intd5,1) );                
        id  p1.p3 = 1/2 * (   1/s2m - 1/s1m - 1/s3m + M^2 );
        id  p1.p6 = 1/2 * (   1/s4m - 1/s1m - p6.p6 );
        id  p3.p6 = 1/2 * ( - 1/s5m + 1/s3m + p6.p6 );
        endif;        
        #call ACCU(D5 2) 

        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

        if ( count(intd5,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p5.p5 = 1/s5m - M^2;
        endif;        
        #call ACCU(D5 3)

        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

        #call ACCU(D5)

* #include matad.info # time
* #include matad.info # print

************************************************************

* do recursion

        #message do recursion

* sort: n3 >= n4,n5,n1
* n3 must be the largest. This is needed for the reduction of
* n6 from <0 to =0. 
        if ( count(intd5,1) );        
        repeat;
                if (count(s3m,1) < count(s4m,1)) 
                multiply replace_(s3m,s4m,s4m,s3m,s1m,s5m,s5m,s1m);
                if (count(s3m,1) < count(s5m,1)) 
                multiply replace_(s3m,s5m,s5m,s3m,s1m,s4m,s4m,s1m);
                if (count(s3m,1) < count(s1m,1)) 
                multiply replace_(s3m,s1m,s1m,s3m,s5m,s4m,s4m,s5m);
        endrepeat;
        endif;        
        .sort

        if ( count(intd5,1) );                
        repeat;
                #call redD5n6p(s1m,s2m,s3m,s4m,s5m,p6)
        endrepeat;
        endif;        
        .sort

        #call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

        #call ACCU(D5)

        #do i=1,1

                #call redD5n6m(s1m,s2m,s3m,s4m,s5m,p6)

* sort: n3 >= n4,n5,n1
* n3 must be the largest. This is needed for the reduction of
* n6 from <0 to =0. 
                if ( count(intd5,1) );        
                repeat;
                        if (count(s3m,1) < count(s4m,1)) 
                        multiply replace_(s3m,s4m,s4m,s3m,s1m,s5m,s5m,s1m);
                        if (count(s3m,1) < count(s5m,1)) 
                        multiply replace_(s3m,s5m,s5m,s3m,s1m,s4m,s4m,s1m);
                        if (count(s3m,1) < count(s1m,1)) 
                        multiply replace_(s3m,s1m,s1m,s3m,s5m,s4m,s4m,s5m);
                endrepeat;
                endif;
                #call ACCU(D5)

#enddo

#call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

* #include expandnomdeno

#call ACCU(D5)

#do i=1,1

        #call redD5n2(s1m,s2m,s3m,s4m,s5m,p6)
        #call ACCU(D5)

#enddo

* sort: n1 >= n3,n4,n5
if ( count(intd5,1) );        
repeat;
        if (count(s1m,1) < count(s3m,1)) 
        multiply replace_(s1m,s3m,s3m,s1m,s4m,s5m,s5m,s4m);
        if (count(s1m,1) < count(s4m,1)) 
        multiply replace_(s1m,s4m,s4m,s1m,s3m,s5m,s5m,s3m);
        if (count(s1m,1) < count(s5m,1)) 
        multiply replace_(s1m,s5m,s5m,s1m,s3m,s4m,s4m,s3m);
endrepeat;
endif;
.sort

#call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

* #include expandnomdeno

#do i=1,1

        #call redD5n1345(s1m,s2m,s3m,s4m,s5m,p6)
        .sort

        id n = num(4-2*ep);
*         repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
        #call ACCU(D5)

#enddo

#call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

id n = num(4-2*ep);
* repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
#call ACCU(D5)

if ( count(intd5,1) );        
if ( (count(s2m,1)==0) && (count(p6,1)==0) );
if (count(s1m,1) < count(s3m,1)) multiply replace_(s1m,s3m,s3m,s1m);
if (count(s1m,1) < count(s4m,1)) multiply replace_(s1m,s4m,s4m,s1m);
if (count(s1m,1) < count(s5m,1)) multiply replace_(s1m,s5m,s5m,s1m);
if (count(s3m,1) < count(s4m,1)) multiply replace_(s3m,s4m,s4m,s3m);
if (count(s3m,1) < count(s5m,1)) multiply replace_(s3m,s5m,s5m,s3m);
if (count(s4m,1) < count(s5m,1)) multiply replace_(s4m,s5m,s5m,s4m);
endif;
if ( (count(s3m,1)==0) && (count(p6,1)==0) );
if (count(s2m,1) < count(s4m,1)) multiply replace_(s2m,s4m,s4m,s2m);
if (count(s2m,1) < count(s5m,1)) multiply replace_(s2m,s5m,s5m,s2m);
if (count(s4m,1) < count(s5m,1)) multiply replace_(s4m,s5m,s5m,s4m);
endif;
* topd5
endif;

#call ACCU(D5)

* transform D5(1,1,1,1,1,0) to D5(1,1,1,1,1,1)
*
* or: compute D5(1,1,1,1,1,0) once and use it as master integral ?

if ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s3m,1)==1) &&
(count(s4m,1)==1) && (count(s5m,1)==1) && (count(p6,1)==0) );
id s1m*s2m*s3m*s4m*s5m = 

M^2*nom(0,-2)*deno(2-8/3,4/3) * (

s1m*s2m*s3m*s4m*s5m/p6.p6 - (

+ s1m^2*s2m*s3m*s5m * (
- 2*deno(0,-2)*p6.p6^-1
)

+ s1m^2*s2m*s4m*s5m * (
+ 2/3*deno(0,-2)*M^-2
)

+ s1m^2*s3m*s4m*s5m * (
- 2/3*deno(0,-2)*M^-2
)
)

);
endif;
.sort

* sort: n5 <= n1,n3,n4
if ( count(intd5,1) );        
repeat;
        if (count(s5m,1) > count(s1m,1)) 
        multiply replace_(s5m,s1m,s1m,s5m,s4m,s3m,s3m,s4m);
        if (count(s5m,1) > count(s3m,1)) 
        multiply replace_(s5m,s3m,s3m,s5m,s1m,s4m,s4m,s1m);
        if (count(s5m,1) > count(s4m,1)) 
        multiply replace_(s5m,s4m,s4m,s5m,s3m,s1m,s1m,s3m);
endrepeat;
endif;
.sort

#call mltadD5(s1m,s2m,s3m,s4m,s5m,p6)

#call ACCU(D5)

************************************************************

* identify simple integrals
if ( count(intd5,1) );        
if ( (count(s2m,1)<=0) );
id 1/s2m=p2.p2+M^2; 
id p2=p1+p3;
id p1=-p1;
id p4=-p4;
id p5=-p5;
multiply replace_(p1,p4,p3,p6,p4,p5,p5,p3,p6,p1,
s1m,s4m,s3m,s6m,s4m,s5m,s5m,s3m);

*** needed because of topBN
id 1/s3m=p3.p3+M^2;

multiply, intbn/intd5;
elseif ( (count(s5m,1)<=0) );
id 1/s5m=p5.p5+M^2; 
id p5=-p5;
multiply replace_(p2,p4,p3,p6,p4,p2,p6,p3,
s2m,s4m,s3m,s6m,s4m,s2m);
multiply, intd4/intd5;
elseif ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s3m,1)==1) && 
(count(s4m,1)==1) && (count(s5m,1)==1) && (count(p6.p6,1)==-1) );
id s1m*s2m*s3m*s4m*s5m/p6.p6 = miD5;
Multiply int0/intd5;
else;
multiply 1/(1-1);
endif;
endif;
*#include expandnomdeno

#message - done
#endprocedure











#procedure mltadD4(s1m,s2m,p3,s4m,p5,s6m)

* discard massless tadpoles
        if ( count(intd4,1) );
        if ( (count(`s1m',1)<=0) && (count(`s6m',1)<=0) ) discard;
        if ( (count(`s4m',1)<=0) && (count(`s6m',1)<=0) ) discard;
        if ( (count(`s2m',1)<=0) && (count(`p3'.`p3',1)>=0) ) discard;
        if ( (count(`s2m',1)<=0) && (count(`p5'.`p5',1)>=0) ) discard;

        if ( (count(`p3'.`p3',1)>=0) && (count(`p5'.`p5',1)>=0) &&
        (count(`s1m',1)<=0)     && (count(`s6m',1)<=0) ) discard;
        endif;
        .sort

#endprocedure

#procedure redD4n35m(s1m,s2m,p3,s4m,p5,s6m)

* reduce n3,n5 from <0 to =0
* (rec. rel. implemented for n3)

        if ( count(intd4,1) );
        if ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m'^2)>0) &&
        (count(p5.p5,1)!=0) );


* for the special case (n3=-2, n5=-1) reverse order of n3 and n5
        if ( (count(`p3'.`p3',1)==2) && (count(`p5'.`p5',1)==1) ) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');


        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n6-1)*(
        - `p3'.`p3'^-1*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        - `p3'.`p3'^-1*`p5'.`p5'*(n6-1)
        + `p3'.`p3'^-1*n4*`s2m'^-1*`s4m'*`s6m'^-1
        - `p3'.`p3'^-1*n4*`s6m'^-1
        - 2*`p3'.`p3'^-1*n5*`s6m'^-1
        + `p3'.`p3'^-1*M^2*(n6-1)
        - `p3'.`p3'^-1*(n6-1)*`s6m'^-1
        + `p3'.`p3'^-1*`s6m'^-1*n
        );

        redefine i "0";

        elseif ( (match(`s1m'^2*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s6m',1)==1) && (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n1-1)*(
        + `p3'.`p3'^-1*`p5'.`p5'*n4*`s1m'^-1*`s4m'
        + 2*`p3'.`p3'^-1*n3*`s1m'^-1
        - `p3'.`p3'^-1*n4*`s1m'^-1*`s2m'^-1*`s4m'
        + `p3'.`p3'^-1*n4*`s1m'^-1
        + 2*`p3'.`p3'^-1*n5*`s1m'^-1
        - 2*`p3'.`p3'^-1*n6*M^2*`s1m'^-1*`s6m'
        + 2*`p3'.`p3'^-1*n6*`s1m'^-1
        + `p3'.`p3'^-1*(n1-1)*`s1m'^-1
        - `p3'.`p3'^-1*(n1-1)*`s2m'^-1
        - 2*`p3'.`p3'^-1*`s1m'^-1*n
        + 2*`p3'.`p3'^-1*`s1m'^-1
        );

        redefine i "0";

        elseif ( (match(`s1m'*`s2m'^2*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s6m',1)==1) && (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n2-1)*(
        + `p3'.`p3'^-1*n1*M^2*`s1m'*`s2m'^-1
        - `p3'.`p3'^-1*n1*`s2m'^-1
        + `p3'.`p3'^-1*n3*`s2m'^-1
        - `p3'.`p3'^-1*n4*`s1m'^-1*`s2m'^-1*`s4m'
        + `p3'.`p3'^-1*n4*`s2m'^-1*`s4m'*`s6m'^-1
        + `p3'.`p3'^-1*n5*`s2m'^-1
        - `p3'.`p3'^-1*n6*M^2*`s2m'^-1*`s6m'
        + `p3'.`p3'^-1*n6*`s2m'^-1
        + `p3'.`p3'^-1*M^2*(n2-1)
        - `p3'.`p3'^-1*(n2-1)*`s1m'^-1
        - 1/2*`p3'.`p3'^-1*`s2m'^-1*n
        + `p3'.`p3'^-1*`s2m'^-1
        );

        redefine i "0";

        elseif ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) && (count(`s6m',1)==1) && 
        (count(`p5'.`p5',1)!=-1) && (count(p5.p5,1)!=0));

*** for the special case (n3=-2, n5=-1) reverse order of n3 and n5
        if ( (count(`p3'.`p3',1)==2) && (count(`p5'.`p5',1)==1) ) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');


        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n5-1)*(
        + `p3'.`p3'^-1*`p5'.`p5'
        - `p3'.`p3'^-1*`p5'.`p5'*n1*M^2*`s1m'
        + `p3'.`p3'^-1*`p5'.`p5'*n1
        - `p3'.`p3'^-1*`p5'.`p5'*n2*M^2*`s2m'
        + `p3'.`p3'^-1*`p5'.`p5'*n2
        + `p3'.`p3'^-1*`p5'.`p5'*n3
        + `p3'.`p3'^-1*`p5'.`p5'*n4*`s1m'^-1*`s4m'
        - `p3'.`p3'^-1*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        + `p3'.`p3'^-1*`p5'.`p5'*n6*M^2*`s6m'
        - `p3'.`p3'^-1*`p5'.`p5'*n6
        - 1/2*`p3'.`p3'^-1*`p5'.`p5'*n
        + `p3'.`p3'^-1*M^2*(n5-1)
        - `p3'.`p3'^-1*(n5-1)*`s6m'^-1
        );

        redefine i "0";

        elseif ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) && (count(`s6m',1)==1) && 
        (count(`p5'.`p5',1)==-1) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)*deno(4-2*n5-n2-n3,-2)*(
        + n2*`s2m'*`s4m'^-1
        - `p3'.`p3'^-1*`p5'.`p5'*n3
        - `p3'.`p3'^-1*n3*M^2
        + `p3'.`p3'^-1*n3*`s6m'^-1
        - `p5'.`p5'*n2*`s2m'
        );

        redefine i "0";

        endif;

        if ( ( (count(`p3',1) < count(`p5',1)) ) && (count(`p3',1)<0) )
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');

* special case (n5=-1):
        if ( (count(`p5'.`p5',1)==1) )
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure

****************************************

* 24.07.98: new realization of "redD4n35m"; not (yet) used 

#procedure redD4n35ma(s1m,s2m,p3,s4m,p5,s6m)

* reduce n3,n5 from <0 to =0
* (rec. rel. implemented for n3)
        if ( count(intd4,1) );
        if ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m'^2)>0) &&
        (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n6-1)*(
        - `p3'.`p3'^-1*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        - `p3'.`p3'^-1*`p5'.`p5'*(n6-1)
        + `p3'.`p3'^-1*n4*`s2m'^-1*`s4m'*`s6m'^-1
        - `p3'.`p3'^-1*n4*`s6m'^-1
        - 2*`p3'.`p3'^-1*n5*`s6m'^-1
        + `p3'.`p3'^-1*M^2*(n6-1)
        - `p3'.`p3'^-1*(n6-1)*`s6m'^-1
        + `p3'.`p3'^-1*`s6m'^-1*n
        );

        elseif ( (match(`s1m'^2*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s6m',1)==1) && (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

*** for this special case (n3=-2, n5=-1) reverse order of n3 and n5
        if ( (count(`p3'.`p3',1)==2) && (count(`p5'.`p5',1)==1) ) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');

        (-1)/(n1-1)*(
        + `p3'.`p3'^-1*`p5'.`p5'*n4*`s1m'^-1*`s4m'
        + 2*`p3'.`p3'^-1*n3*`s1m'^-1
        - `p3'.`p3'^-1*n4*`s1m'^-1*`s2m'^-1*`s4m'
        + `p3'.`p3'^-1*n4*`s1m'^-1
        + 2*`p3'.`p3'^-1*n5*`s1m'^-1
        - 2*`p3'.`p3'^-1*n6*M^2*`s1m'^-1*`s6m'
        + 2*`p3'.`p3'^-1*n6*`s1m'^-1
        + `p3'.`p3'^-1*(n1-1)*`s1m'^-1
        - `p3'.`p3'^-1*(n1-1)*`s2m'^-1
        - 2*`p3'.`p3'^-1*`s1m'^-1*n
        + 2*`p3'.`p3'^-1*`s1m'^-1
        );

        elseif ( (match(`s1m'*`s2m'^2*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s6m',1)==1) && (count(p5.p5,1)!=0) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n2-1)*(
        + `p3'.`p3'^-1*n1*M^2*`s1m'*`s2m'^-1
        - `p3'.`p3'^-1*n1*`s2m'^-1
        + `p3'.`p3'^-1*n3*`s2m'^-1
        - `p3'.`p3'^-1*n4*`s1m'^-1*`s2m'^-1*`s4m'
        + `p3'.`p3'^-1*n4*`s2m'^-1*`s4m'*`s6m'^-1
        + `p3'.`p3'^-1*n5*`s2m'^-1
        - `p3'.`p3'^-1*n6*M^2*`s2m'^-1*`s6m'
        + `p3'.`p3'^-1*n6*`s2m'^-1
        + `p3'.`p3'^-1*M^2*(n2-1)
        - `p3'.`p3'^-1*(n2-1)*`s1m'^-1
        - 1/2*`p3'.`p3'^-1*`s2m'^-1*n
        + `p3'.`p3'^-1*`s2m'^-1
        );

        endif;

        if (count(`p3',1) < count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure

#procedure redD4n35md(s1m,s2m,p3,s4m,p5,s6m)
        if ( count(intd4,1) );
        if ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) && (count(`s6m',1)==1) && 
        (count(`p5'.`p5',1)!=-1) && (count(`p5'.`p5',1)!=0));

* for this special case (n3=-2, n5=-1) reverse order of n3 and n5
        if ( (count(`p3'.`p3',1)==2) && (count(`p5'.`p5',1)==1) ) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/(n5-1)*(
        + `p3'.`p3'^-1*`p5'.`p5'
        - `p3'.`p3'^-1*`p5'.`p5'*n1*M^2*`s1m'
        + `p3'.`p3'^-1*`p5'.`p5'*n1
        - `p3'.`p3'^-1*`p5'.`p5'*n2*M^2*`s2m'
        + `p3'.`p3'^-1*`p5'.`p5'*n2
        + `p3'.`p3'^-1*`p5'.`p5'*n3
        + `p3'.`p3'^-1*`p5'.`p5'*n4*`s1m'^-1*`s4m'
        - `p3'.`p3'^-1*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        + `p3'.`p3'^-1*`p5'.`p5'*n6*M^2*`s6m'
        - `p3'.`p3'^-1*`p5'.`p5'*n6
        - 1/2*`p3'.`p3'^-1*`p5'.`p5'*n
        + `p3'.`p3'^-1*M^2*(n5-1)
        - `p3'.`p3'^-1*(n5-1)*`s6m'^-1
        );

        endif;

        if (count(`p3',1) < count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure


#procedure redD4n35me(s1m,s2m,p3,s4m,p5,s6m)
        if ( count(intd4,1) );
        if ( (match(`s1m'*`s2m'*`p3'.`p3'*`s4m'          *`s6m')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) && (count(`s6m',1)==1) && 
        (count(`p5'.`p5',1)==-1) );

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)*deno(4-2*n5-n2-n3,-2)*(
        + n2*`s2m'*`s4m'^-1
        - `p3'.`p3'^-1*`p5'.`p5'*n3
        - `p3'.`p3'^-1*n3*M^2
        + `p3'.`p3'^-1*n3*`s6m'^-1
        - `p5'.`p5'*n2*`s2m'
        );

        endif;

        if (count(`p3',1) < count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure

****************************************

#procedure redD4n35(s1m,s2m,p3,s4m,p5,s6m)

* reduce n3,n5 from >1 to =1
* (rec. rel. implemented for n3)

        if ( count(intd4,1) );
        if (match(`s1m'*`s2m'/`p3'.`p3'^2*`s4m'/`p5'.`p5'*`s6m')>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/M^2/(n3-1)*(
        - (n3-1)*`s6m'^-1
        + `p3'.`p3'*`p5'.`p5'*n2*`s2m'
        - `p3'.`p3'*n2*`s2m'*`s4m'^-1
        + `p3'.`p3'*n2
        + 2*`p3'.`p3'*n5
        + `p3'.`p3'*(n3-1)
        - `p3'.`p3'*n
        + `p5'.`p5'*(n3-1)
        );

        redefine i "0";

        endif;

        if (count(`p3',1) > count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure


#procedure redD4n2(s1m,s2m,p3,s4m,p5,s6m)

* reduce n2 from >1 to =1
        if ( count(intd4,1) );
        if (match(`s1m'*`s2m'^2/`p3'.`p3'*`s4m'/`p5'.`p5'*`s6m')>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/M^2/(n2-1)*(
        - 1/2*n3*`s2m'^-1
        - 1/2*n5*`s2m'^-1
        - (n2-1)*`s2m'^-1
        + 1/2*`s2m'^-1*n
        + 1/2*`p3'.`p3'^-1*n3*`s1m'^-1*`s2m'^-1
        - 1/2*`p3'.`p3'^-1*n3*`s2m'^-2
        - 1/2*`p5'.`p5'^-1*n5*`s2m'^-2
        + 1/2*`p5'.`p5'^-1*n5*`s2m'^-1*`s4m'^-1
        );

        redefine i "0";

        endif;

        if (count(`p3',1) > count(`p5',1)) 
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        endif;
#endprocedure


#procedure redD4n14(s1m,s2m,p3,s4m,p5,s6m)

* reduce n1,n4 from >1 to =1
* (rec. rel. implemented for n1)
        if ( count(intd4,1) );
        if (match(`s1m'^2*`s2m'/`p3'.`p3'*`s4m'/`p5'.`p5'*`s6m')>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/M^2/(n1-1)*(
        + n2*`s1m'^-1*`s2m'*`s4m'^-1
        - n2*`s1m'^-1
        - n4*`s1m'^-1*`s2m'^-1*`s4m'
        + n4*`s1m'^-1
        - (n1-1)*`s2m'^-1
        + (n1-1)*`s4m'^-1
        - (n1-1)*`s6m'^-1
        + `p3'.`p3'*(n1-1)
        - `p5'.`p5'*n2*`s1m'^-1*`s2m'
        + `p5'.`p5'*n4*`s1m'^-1*`s4m'
        );

        redefine i "0";

        endif;

***if (count(`s4m',1) > count(`s1m',1)) 
***  multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
*** 11Mar04: Changed to (see email from M. Faisst)
        if (count(`s4m',1) > count(`s1m',1)) ;
        multiply replace_(`s1m',`s4m',`s4m',`s1m',`p3',`p5',`p5',`p3');
        redefine i "0" ;
        endif;
        endif;
#endprocedure


#procedure redD4n6(s1m,s2m,p3,s4m,p5,s6m)

* reduce n6 from >1 to =1
        if ( count(intd4,1) );
        if (match(`s1m'*`s2m'/`p3'.`p3'*`s4m'/`p5'.`p5'*`s6m'^2)>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?*`s4m'^n4?/`p5'.`p5'^n5?*`s6m'^n6? =
        `s1m'^n1 *`s2m'^n2 /`p3'.`p3'^n3 *`s4m'^n4 /`p5'.`p5'^n5 *`s6m'^n6 * 

        (-1)/M^2/(n6-1)*(
        + 1/2*n1*`s1m'*`s2m'^-1*`s6m'^-1
        - 1/2*n1*`s6m'^-1
        - n3*`s6m'^-1
        + 1/2*n4*`s2m'^-1*`s4m'*`s6m'^-1
        - 1/2*n4*`s6m'^-1
        - n5*`s6m'^-1
        - (n6-1)*`s6m'^-1
        + `s6m'^-1*n
        - 1/2*`p3'.`p3'*n1*`s1m'*`s6m'^-1
        - 1/2*`p5'.`p5'*n4*`s4m'*`s6m'^-1
        );

        redefine i "0";

        endif;
        endif;
#endprocedure



#procedure topd4
*
* this is topd4
*
* 20.Jun.02: splitting of diad4 commented
*
        #message this is topd4

************************************************************
*
* Jul. 24th 1998:
* - split expression before the repeat-endrepeat-loop of "redD4n35m"
* - split furthermore the procedure "redD4n35m" into smaller parts 
*   (not extensively tested)
*
************************************************************
        
************************************************************

* treat the scalar products

        #message numerator
        if ( count(intd4,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p6.p6 = 1/s6m - M^2;
        endif;        
        #call ACCU(D4 0)

        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

        if ( count(intd4,1) );                
        id p2 = p1+p3;
        id p4 = p1+p6;
        id p5 = p6-p3;
        endif;        
        #call ACCU(D4 1)

        if ( count(intd4,1) );                
        id  p1.p3 = 1/2 * (   1/s2m - 1/s1m - p3.p3 );
        id  p1.p6 = 1/2 * (   1/s4m - 1/s1m - 1/s6m + M^2 );
        id  p3.p6 = 1/2 * ( - p5.p5 + p3.p3 + 1/s6m - M^2 );
        endif;        
        #call ACCU(D4 2)

        if ( count(intd4,1) );                
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p4.p4 = 1/s4m - M^2;
        id p6.p6 = 1/s6m - M^2;
        endif;        
        #call ACCU(D4 3)

        #call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

        #call ACCU(D4)

* #include matad.info # time
* #include matad.info # print

************************************************************
*
* do recursion
*

        #message do recursion

*
* massless indices <0
*
        #do i=1,1

                #call redD4n35m(s1m,s2m,p3,s4m,p5,s6m)
                #call ACCU(D4n35m)

#enddo

#call mltadD4(s1m,s2m,p3,s4m,p5,s6m)
* #include expandnomdeno

**** As "redD4n35m" has a complicated structure the terms which are
**** already reduced have to be split from the rest.

***G diad4h1 = diad4;
***.sort
***drop diad4;
***.sort
***hide diad4h1;
***.sort

***G diad4h2 = diad4h1;
***if ( (match(s1m*s2m*s4m*s6m)>0) && 
***   ( ( (count(p3.p3,1)>0) && (count(p5.p5,1)!=0) ) ||
***     ( (count(p5.p5,1)>0) && (count(p3.p3,1)!=0) ) ) ) discard;
***.sort
***hide diad4h2;
***.sort

***G diad4h3 = diad4h1 - diad4h2;
***.sort

#do i=1,1

        #call redD4n35m(s1m,s2m,p3,s4m,p5,s6m)
        #call ACCU(D4n35m)

#enddo

#call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

* AFP
* #include expandnomdeno

***G diad4 = diad4h2 + diad4h3;
***.sort

***unhide diad4h1,diad4h2,diad4h3;
***.sort
***drop diad4h1,diad4h2,diad4h3;
***.sort

*
* massless indices >0
*

#do i=1,1

        #call redD4n35(s1m,s2m,p3,s4m,p5,s6m)
        #call redD4n2(s1m,s2m,p3,s4m,p5,s6m)
        #call ACCU(D4n235)

#enddo

#call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

* #include expandnomdeno

#do i=1,1

        #call redD4n14(s1m,s2m,p3,s4m,p5,s6m)
        .sort

#enddo

#call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

#call ACCU(D4)

#do i=1,1

        #call redD4n35(s1m,s2m,p3,s4m,p5,s6m)
        #call redD4n2(s1m,s2m,p3,s4m,p5,s6m)
        #call ACCU(D4n235)

#enddo

#call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

* #include expandnomdeno

#do i=1,1

        #call redD4n6(s1m,s2m,p3,s4m,p5,s6m)
        .sort

#enddo

#call mltadD4(s1m,s2m,p3,s4m,p5,s6m)

#call ACCU(D4)

if ( count(intd4,1) );        
if (count(s4m,1) < count(s1m,1)) 
multiply replace_(s1m,s4m,s4m,s1m,p3,p5,p5,p3);

if ( (count(s1m,1)>0) && (count(s4m,1)>0) );
if (count(p3,1) < count(p5,1)) 
multiply replace_(s1m,s4m,s4m,s1m,p3,p5,p5,p3);
endif;

if (count(p5,1)==0) multiply replace_(s1m,s4m,s4m,s1m,p3,p5,p5,p3);
endif;
.sort

#call ACCU(D4)

************************************************************

* identify simple integrals

if ( count(intd4,1) );        
if ( (count(s6m,1)<=0) && (count(s2m,1)>0) );
id 1/s6m=p6.p6+M^2; 
id p6=p4-p1;
id p5=-p5;
multiply replace_(p1,p5,p2,p4,p3,p2,p4,p6,p5,p1,
s1m,s5m,s2m,s4m,s4m,s6m);
multiply, intbm/intd4;
elseif ( (count(s2m,1)<=0) );
id 1/s2m=p2.p2+M^2;
id p2=p1+p3; 
id p3=-p3;
id p4=-p4;
id p6=-p6;
multiply replace_(p3,p4,p4,p3,p6,p2,
s4m,s3m,s6m,s2m);
multiply, intdm/intd4;
elseif ( (count(s1m,1)<=0) );
id 1/s1m=p1.p1+M^2;
id p1=p4-p6; 
id p5=-p5;
multiply replace_(p2,p4,p3,p5,p4,p6,p5,p2,p6,p3,
s2m,s4m,s4m,s6m,s6m,s3m);

* needed because of topBN1
id 1/s3m=p3.p3+M^2;
id 1/s4m=p4.p4+M^2;
id 1/s6m=p6.p6+M^2;

multiply, intbn1/intd4;
elseif ( (count(p3,1)==0) );
id p5=-p5;
multiply replace_(p1,p3,p2,p4,p4,p2,p6,p1,
s1m,s3m,s2m,s4m,s4m,s2m,s6m,s1m);
multiply, inte4/intd4;
elseif ( (count(s1m,1)==1) && (count(s2m,1)==1)    && (count(p3.p3,1)==-1) && 
(count(s4m,1)==1) && (count(p5.p5,1)==-1) && (count(s6m,1)==1) );
id s1m*s2m/p3.p3*s4m/p5.p5*s6m = miD4;
Multiply int0/intd4;
else;
multiply 1/(1-1);
endif;
endif;
#call ACCU(D4)

#message - done
#endprocedure






************************************************************

#procedure mltadDM(s1m,s2m,s3m,p4,p5,p6)

* discard massless tadpoles
*
* two massless indices >=0 => D_M = 0
        if ( count(intdm,1) );                
        if ( (count(`p6'.`p6',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;
        if ( (count(`p6'.`p6',1)>=0) && (count(`p5'.`p5',1)>=0) ) discard;
        if ( (count(`p5'.`p5',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;

* two massive indices >=0 => D_M = 0

        if ( (count(`s1m',1)<=0) && (count(`s2m',1)<=0) ) discard;
        if ( (count(`s1m',1)<=0) && (count(`s3m',1)<=0) ) discard;
        if ( (count(`s2m',1)<=0) && (count(`s3m',1)<=0) ) discard;
        endif;        
        .sort

#endprocedure

#procedure redDMn4m(s1m,s2m,s3m,p4,p5,p6)

* reduce massless indices from <0 to =0

* sort: n5,n6 >= n4
        if ( count(intdm,1) );                
        if (count(`p5'.`p5',1) > count(`p4'.`p4',1)) 
        multiply replace_(`p4',`p5',`p5',`p4',`s1m',`s3m',`s3m',`s1m');
        if (count(`p4'.`p4',1) < count(`p6'.`p6',1)) 
        multiply replace_(`p6',`p4',`p4',`p6',`s3m',`s2m',`s2m',`s3m');


*** n1=n2=1:
        if ( (match(`s1m'*`s2m'*`s3m'*`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0) &&
        (count(`s1m',1)==1) && (count(`s2m',1)==1) );

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)*deno(4-n4-3/2*n1-3/2*n2,-2)*(
        + 3/2*n1*M^2*`s1m'
        - 1/2*n1*`s1m'*`s2m'^-1
        + 1/2*n1*`s1m'*`s3m'^-1
        + 3/2*n2*M^2*`s2m'
        - 1/2*n2*`s1m'^-1*`s2m'
        + 1/2*n2*`s2m'*`s3m'^-1
        + 1/2*`p4'.`p4'^-1*`p5'.`p5'*n4
        + 1/2*`p4'.`p4'^-1*`p6'.`p6'*n4
        + `p4'.`p4'^-1*n4*M^2
        - 1/2*`p4'.`p4'^-1*n4*`s1m'^-1
        - 1/2*`p4'.`p4'^-1*n4*`s2m'^-1
        );

        redefine i "0";

        endif;

        if (match(`s1m'*`s2m'^2*`s3m'*`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)/(n2-1)*(
        - `p4'.`p4'^-1*`p5'.`p5'*n3*`s2m'^-1*`s3m'
        - `p4'.`p4'^-1*`p5'.`p5'*(n2-1)
        + `p4'.`p4'^-1*`p6'.`p6'*n3*`s2m'^-1*`s3m'
        - `p4'.`p4'^-1*n1*M^2*`s1m'*`s2m'^-1
        + `p4'.`p4'^-1*n1*`s2m'^-1
        + `p4'.`p4'^-1*n4*`s2m'^-1
        - `p4'.`p4'^-1*n5*`s2m'^-1
        + `p4'.`p4'^-1*n6*`s2m'^-1
        - 1/2*`p4'.`p4'^-1*`s2m'^-1*n
        + `p4'.`p4'^-1*`s2m'^-1
        );

        redefine i "0";

        endif;

        if ( (count(`p5'.`p5',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;
        if ( (count(`p6'.`p6',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;


        if ( (match(`s1m'^2*`s2m'*`s3m'*`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0) &&
        (count(`s2m',1)==1) );

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)/(n1-1)*(
        + `p4'.`p4'^-1*`p5'.`p5'*n3*`s1m'^-1*`s3m'
        - `p4'.`p4'^-1*`p6'.`p6'*n3*`s1m'^-1*`s3m'
        - `p4'.`p4'^-1*`p6'.`p6'*(n1-1)
        - `p4'.`p4'^-1*n2*M^2*`s1m'^-1*`s2m'
        + `p4'.`p4'^-1*n2*`s1m'^-1
        + `p4'.`p4'^-1*n4*`s1m'^-1
        + `p4'.`p4'^-1*n5*`s1m'^-1
        - `p4'.`p4'^-1*n6*`s1m'^-1
        - 1/2*`p4'.`p4'^-1*`s1m'^-1*n
        + `p4'.`p4'^-1*`s1m'^-1
        );

        redefine i "0";

        endif;

        if ( (count(`p5'.`p5',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;
        if ( (count(`p6'.`p6',1)>=0) && (count(`p4'.`p4',1)>=0) ) discard;
        endif;
#endprocedure

#procedure redDMn456(s1m,s2m,s3m,p4,p5,p6)

* reduce massless indices from >1 to =1

* sort: n5 >= n4 >= n6
        if ( count(intdm,1) );                
        if (count(`p5'.`p5',1) > count(`p4'.`p4',1)) 
        multiply replace_(`p4',`p5',`p5',`p4',`s1m',`s3m',`s3m',`s1m');
        if (count(`p5'.`p5',1) > count(`p6'.`p6',1)) 
        multiply replace_(`p6',`p5',`p5',`p6',`s1m',`s2m',`s2m',`s1m');
        if (count(`p4'.`p4',1) > count(`p6'.`p6',1)) 
        multiply replace_(`p6',`p4',`p4',`p6',`s3m',`s2m',`s2m',`s3m');

        if (match(`s1m'*`s2m'*`s3m'/`p4'.`p4'/`p5'.`p5'^2/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 
        (-1)/M^2/(n5-1)*(
        - 1/2*(n5-1)*`s2m'^-1
        - 1/2*(n5-1)*`s3m'^-1
        - 3/2*`p4'.`p4'*`p5'.`p5'*n2*`s2m'
        + 1/2*`p4'.`p4'*(n5-1)
        + 3/2*`p5'.`p5'^2*n2*`s2m'
        + 3/2*`p5'.`p5'^2*n3*`s3m'
        - 3/2*`p5'.`p5'*`p6'.`p6'*n3*`s3m'
        + 1/2*`p5'.`p5'*n2*`s1m'^-1*`s2m'
        - 1/2*`p5'.`p5'*n2*`s2m'*`s3m'^-1
        + 1/2*`p5'.`p5'*n3*`s1m'^-1*`s3m'
        - 1/2*`p5'.`p5'*n3*`s2m'^-1*`s3m'
        + 2*`p5'.`p5'*(n5-1)
        - 1/2*`p5'.`p5'*n
        + 1/2*`p6'.`p6'*(n5-1)
        );

        redefine i "0";

        endif;
        endif;
#endprocedure

#procedure redDMn123(s1m,s2m,s3m,p4,p5,p6)

* reduce massive indices from >1 to =1

* sort: n1 >= n2 >= n3
        if ( count(intdm,1) );                
        if (count(`s1m',1) < count(`s2m',1)) 
        multiply replace_(`p6',`p5',`p5',`p6',`s1m',`s2m',`s2m',`s1m');
        if (count(`s1m',1) < count(`s3m',1)) 
        multiply replace_(`p4',`p5',`p5',`p4',`s1m',`s3m',`s3m',`s1m');
        if (count(`s2m',1) < count(`s3m',1)) 
        multiply replace_(`p6',`p4',`p4',`p6',`s3m',`s2m',`s2m',`s3m');

        if (match(`s1m'^2*`s2m'*`s3m'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?/`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 
        (-1)/M^2/(n1-1)*(
        - n4*`s1m'^-1
        + n5*`s1m'^-1
        - n6*`s1m'^-1
        - (n1-1)*`s1m'^-1
        + 1/2*`s1m'^-1*n
        - `p4'.`p4'*n2*`s1m'^-1*`s2m'
        + `p5'.`p5'*n2*`s1m'^-1*`s2m'
        + `p5'.`p5'*n3*`s1m'^-1*`s3m'
        - `p6'.`p6'*n3*`s1m'^-1*`s3m'
        );

        redefine i "0";

        endif;
        endif;
#endprocedure




#procedure topdm

*
* this is topdm
*
        #message this is topdm


************************************************************

* treat the scalar products

        #message numerator

        if ( count(intdm,1) );        
        id p2=p1+p3;
        id p6=p3+p5;
        endif;        
        #call ACCU(DM 1)

        if ( count(intdm,1) );        
        id  p1.p3 = 1/2 * (   p2.p2 - p1.p1 - p3.p3 );
        id  p1.p4 = 1/2 * ( - p6.p6 + p1.p1 + p4.p4 );
        endif;        
        #call ACCU(DM 2)

        if ( count(intdm,1) );                
        id  p1.p5 = 1/2 * (   p3.p3 + p4.p4 - p2.p2 - p6.p6 );
        id  p3.p4 = 1/2 * (   p2.p2 + p6.p6 - p1.p1 - p5.p5 );
        endif;        
        #call ACCU(DM 3)

        if ( count(intdm,1) );                
        id  p3.p5 = 1/2 * (   p6.p6 - p3.p3 - p5.p5 );
        id  p4.p5 = 1/2 * ( - p2.p2 + p4.p4 + p5.p5 );
        endif;        
        #call ACCU(DM 4)

        if ( count(intdm,1) );                
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        endif;        
        #call ACCU(DM 5)

        #call mltadDM(s1m,s2m,s3m,p4,p5,p6)

        #call ACCU(DM)

* #include matad.info # time
* #include matad.info # print

************************************************************

* do recursion

        #message do recursion

* sort: n5 >= n4 >= n6

        if ( count(intdm,1) );                
        if (count(p5.p5,1) > count(p4.p4,1)) 
        multiply replace_(p4,p5,p5,p4,s1m,s3m,s3m,s1m);
        if (count(p5.p5,1) > count(p6.p6,1)) 
        multiply replace_(p6,p5,p5,p6,s1m,s2m,s2m,s1m);
        if (count(p4.p4,1) > count(p6.p6,1)) 
        multiply replace_(p6,p4,p4,p6,s3m,s2m,s2m,s3m);
        endif;        
        .sort

************************************************************

* massless index <0

        #do i=1,1
                #call redDMn4m(s1m,s2m,s3m,p4,p5,p6)
                .sort
        #enddo

        #call mltadDM(s1m,s2m,s3m,p4,p5,p6)
* #include expandnomdeno

        id n = num(4-2*ep);
*         id acc(x1?)*acc(x2?) = acc(x1*x2);
        #call ACCU(DM)

* massless index >1

        #do i=1,1
                #call redDMn456(s1m,s2m,s3m,p4,p5,p6)
                .sort
        #enddo

        #call mltadDM(s1m,s2m,s3m,p4,p5,p6)

        id n = num(4-2*ep);
*         id acc(x1?)*acc(x2?) = acc(x1*x2);
        #call ACCU(DM)

        #do i=1,1
                #call redDMn123(s1m,s2m,s3m,p4,p5,p6)
                .sort
        #enddo

        #call mltadDM(s1m,s2m,s3m,p4,p5,p6)

        id n = num(4-2*ep);
*         id acc(x1?)*acc(x2?) = acc(x1*x2);
        #call ACCU(DM)

************************************************************

* identify simple integrals

* sort: n5 >= n4 >= n6
* (this order is needed for the identification of the simple integrals)

        if ( count(intdm,1) );        
        if (count(p5.p5,1) > count(p4.p4,1)) 
        multiply replace_(p4,p5,p5,p4,s1m,s3m,s3m,s1m);
        if (count(p5.p5,1) > count(p6.p6,1)) 
        multiply replace_(p6,p5,p5,p6,s1m,s2m,s2m,s1m);
        if (count(p4.p4,1) > count(p6.p6,1)) 
        multiply replace_(p6,p4,p4,p6,s3m,s2m,s2m,s3m);
        endif;
        .sort

* if one of the massive indices is absent
* sort: n1 <= n2 <= n3

        if ( count(intdm,1) );        
        if (match(s1m*s2m*s3m)==0);
        if (count(s1m,1) > count(s2m,1)) 
        multiply replace_(p6,p5,p5,p6,s1m,s2m,s2m,s1m);
        if (count(s1m,1) > count(s3m,1)) 
        multiply replace_(p4,p5,p5,p4,s1m,s3m,s3m,s1m);
        if (count(s2m,1) > count(s3m,1)) 
        multiply replace_(p6,p4,p4,p6,s3m,s2m,s2m,s3m);
        endif;
        endif;        
        .sort

        #call ACCU(DM)

        if ( count(intdm,1) );        
        if ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2; 
        id p1=p2-p3;
        id p5=-p5;
        multiply replace_(p2,p6,p3,p4,p4,p3,p5,p1,p6,p5,s2m,s6m,s3m,s4m);
        multiply, intbn2/intdm;
        elseif ( (count(s1m,1) >0)   && (count(s2m,1) >0)   && (count(s3m,1) >0) && 
        (count(p4.p4,1) <0) && (count(p5.p5,1) <0) && (count(p6.p6,1)==0) ); 

* use symmetry for E_3:

        if (count(s1m,1) > count(s3m,1))
        multiply replace_(p1,p3,p3,p1,s1m,s3m,s3m,s1m);
        multiply, inte3/intdm;
        elseif ( (count(s1m,1)==1)   && (count(s2m,1)==1)   && (count(s3m,1)==1) && 
        (count(p4.p4,1)==-1) && (count(p5.p5,1)==-1) && (count(p6.p6,1)==-1));
        id s1m*s2m*s3m/p4.p4/p5.p5/p6.p6 = 2*z3/ep + DM + ep*DMep;
        else;
        multiply 1/(1-1);
        endif;
        endif;
        #call ACCU(DM)

        #message - done
        
#endprocedure









************************************************************

#procedure mltadDN(s1m,s2m,p3,p4,p5,p6)

* discard massless tadpoles
        if ( count(intdn,1) );
        if ( (count(`s1m',1)<=0) && (count(`s2m',1)<=0) ) discard;
        if ( (count(`s1m',1)<=0) && 
        ( (count(`p3'.`p3',1)>=0) || (count(`p4'.`p4',1)>=0) ||
        (count(`p5'.`p5',1)>=0) || (count(`p6'.`p6',1)>=0) ) ) discard;
        if ( (count(`s2m',1)<=0) && 
        ( (count(`p3'.`p3',1)>=0) || (count(`p4'.`p4',1)>=0) ||
        (count(`p5'.`p5',1)>=0) || (count(`p6'.`p6',1)>=0) ) ) discard;
        endif;        
        .sort

#endprocedure


#procedure redDNn6(s1m,s2m,p3,p4,p5,p6)

* reduce massless indices from >1 to =1
        if ( count(intdn,1) );
        if (match(`s1m'*`s2m'/`p3'.`p3'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

* sort: n6 >=  n3,n4,n5

        if (count(`p6'.`p6',1) > count(`p5'.`p5',1)) 
        multiply replace_(`p6',`p5',`p5',`p6',`s1m',`s2m',`s2m',`s1m');
        if (count(`p6'.`p6',1) > count(`p4'.`p4',1)) 
        multiply replace_(`p6',`p4',`p4',`p6',`p3',`p5',`p5',`p3');
        if (count(`p6'.`p6',1) > count(`p3'.`p3',1)) 
        multiply replace_(`p6',`p3',`p3',`p6',`p4',`p5',`p5',`p4');
        endif;

        if (match(`s1m'*`s2m'/`p3'.`p3'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6'^2)>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?
        /`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2/`p3'.`p3'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)/M^2/(n6-1)*(
        - 1/2*(n6-1)*`s1m'^-1
        - 1/2*(n6-1)*`s2m'^-1
        + 1/2*`p3'.`p3'*`p6'.`p6'*n2*`s2m'
        + 1/2*`p3'.`p3'*(n6-1)
        + 1/2*`p4'.`p4'*`p6'.`p6'*n1*`s1m'
        + 1/2*`p4'.`p4'*(n6-1)
        - 1/2*`p5'.`p5'*`p6'.`p6'*n1*`s1m'
        - 1/2*`p5'.`p5'*`p6'.`p6'*n2*`s2m'
        + 1/2*`p6'.`p6'*n3
        + 1/2*`p6'.`p6'*n4
        - 1/2*`p6'.`p6'*n5
        + 1/2*`p6'.`p6'*(n6-1)
        - 1/4*`p6'.`p6'*n
        );

        redefine i "0";
        redefine ii "0";
        endif;
        endif;
#endprocedure

#procedure redDNn1(s1m,s2m,p3,p4,p5,p6)

* reduce massive indices from >1 to =1
        if ( count(intdn,1) );
        if (match(`s1m'*`s2m'/`p3'.`p3'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

* sort: n1 >= n2

        if (count(`s1m',1) < count(`s2m',1)) 
        multiply replace_(`p6',`p5',`p5',`p6',`s1m',`s2m',`s2m',`s1m');
        endif;

        if (match(`s1m'^2*`s2m'/`p3'.`p3'/`p4'.`p4'/`p5'.`p5'/`p6'.`p6')>0);

        id `s1m'^n1?*`s2m'^n2?/`p3'.`p3'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5?
        /`p6'.`p6'^n6? =
        `s1m'^n1*`s2m'^n2/`p3'.`p3'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5/`p6'.`p6'^n6 * 

        (-1)/M^2/(n1-1)*(
        + 1/2*n3*`s1m'^-1
        - 3/2*n4*`s1m'^-1
        - 1/2*n5*`s1m'^-1
        - 1/2*n6*`s1m'^-1
        - (n1-1)*`s1m'^-1
        + 3/4*`s1m'^-1*n
        + 1/2*`p3'.`p3'*`p6'.`p6'^-1*n6*`s1m'^-1
        + 1/2*`p3'.`p3'*n2*`s1m'^-1*`s2m'
        - 1/2*`p4'.`p4'*`p6'.`p6'^-1*n6*`s1m'^-1
        - 1/2*`p4'.`p4'*(n1-1)
        - 1/2*`p5'.`p5'*n2*`s1m'^-1*`s2m'
        + 1/2*`p5'.`p5'*(n1-1)
        - 1/2*`p6'.`p6'^-1*n6*`s1m'^-2
        + 1/2*`p6'.`p6'^-1*n6*`s1m'^-1*`s2m'^-1
        );

        redefine i "0";
        redefine ii "0" ;
        endif;
        endif;
#endprocedure



#procedure topdn
*
* this is topdn
*
        #message this is topdn


************************************************************

* treat the scalar products     CHECK !!!!!!!!!!

        #message numerator
        if ( count(intdn,1) );
        id p3=p5-p2;
        id p4=p2+p6;
        endif;        
        #call ACCU(DN 1)

        if ( count(intdn,1) );        
        id  p1.p2 = 1/2 * ( p4.p4 + p3.p3 - p5.p5 - p6.p6 );
        id  p1.p5 = 1/2 * ( p4.p4 - p5.p5 - p1.p1 );
        endif;        
        #call ACCU(DN 2)

        if ( count(intdn,1) );        
        id  p1.p6 = 1/2 * (-p3.p3 + p6.p6 + p1.p1 );
        id  p2.p5 = 1/2 * (-p3.p3 + p5.p5 + p2.p2 );
        endif;        
        #call ACCU(DN 3)

        if ( count(intdn,1) );        
        id  p2.p6 = 1/2 * ( p4.p4 - p6.p6 - p2.p2 );
        id  p5.p6 = 1/2 * ( p3.p3 + p4.p4 - p2.p2 - p1.p1);
        endif;        
        #call ACCU(DN 4)

        if ( count(intdn,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        endif;        
        .sort

        #call mltadDN(s1m,s2m,p3,p4,p5,p6)

        #call ACCU(DN)

* #include matad.info # time
* #include matad.info # print

************************************************************

* do recursion

        #message do recursion

        #define ii "0"
        #do i=1,1
                #call redDNn6(s1m,s2m,p3,p4,p5,p6)
                .sort
        #enddo
        #undefine ii

        #call mltadDN(s1m,s2m,p3,p4,p5,p6)

* #include expandnomdeno

        id n=num(4-2*ep);
*         id acc(x1?)*acc(x2?)=acc(x1*x2);
        #call ACCU(DN)

        #do i=1,1

                #do ii=1,1
                        #call redDNn1(s1m,s2m,p3,p4,p5,p6)
                        .sort
                #enddo 

                #call mltadDN(s1m,s2m,p3,p4,p5,p6)

                #do ii=1,1
                        #call redDNn6(s1m,s2m,p3,p4,p5,p6)
                        .sort
                #enddo

                id n=num(4-2*ep);
*                 id acc(x1?)*acc(x2?)=acc(x1*x2);
                #call ACCU(DN)

#enddo

#call mltadDN(s1m,s2m,p3,p4,p5,p6)

id n=num(4-2*ep);
* id acc(x1?)*acc(x2?)=acc(x1*x2);
#call ACCU(DN)

************************************************************

* identify simple integrals

if ( count(intdn,1) );
if (count(s2m,1)<=0) multiply replace_(p6,p5,p5,p6,s1m,s2m,s2m,s1m);
if (count(p4.p4,1)>=0) multiply replace_(p4,p3,p3,p4,s1m,s2m,s2m,s1m);
if (count(p5.p5,1)>=0) multiply replace_(p5,p3,p3,p5,p4,p6,p6,p4);
if (count(p6.p6,1)>=0) multiply replace_(p6,p3,p3,p6,p4,p5,p5,p4);
endif;
.sort

if ( count(intdn,1) );
if ( (count(s1m,1)<=0) );
id 1/s1m=p1.p1+M^2; 
id p1=p4-p5;
id p3=-p3;
id p5=-p5;
multiply replace_(p2,p6,p4,p2,p5,p4,p6,p1,
s2m,s6m);
multiply, intm3/intdn;
elseif ( (count(p3.p3,1)>=0) );
id p3=p5-p2;
id p6=-p6;
multiply replace_(p1,p3,p4,p5,p2,p4,p5,p2,p6,p1,
s1m,s3m,s2m,s4m);
multiply, intbn1/intdn;
elseif ( (count(s1m,1)==1)   && (count(s2m,1)==1)   && (count(p3.p3,1)==-1) && 
(count(p4.p4,1)==-1) && (count(p5.p5,1)==-1) && (count(p6.p6,1)==-1));
id s1m*s2m/p3.p3/p4.p4/p5.p5/p6.p6 = miDN;
else;
multiply 1/(1-1);
endif;
endif;
#call ACCU(DN)

#message - done

#endprocedure        









#procedure tope4

*
* this is tope4
*
        #message this is tope4

* As this topology results form other ones by shrinking one line
* in principle no scalar products may appear. However, one has to
* take care, that, e.g., in the topology D5 the index of line 2 is 
* really reduced to 0 and not to negative values. Otherwise unwanted 
* scalar products appear.
* Remark: It can be avoided to reduce D5 to E4 ...
        if ( count(inte4,1) );
        if ( count(p1.p5,1,p1.p4,1,p3.p4,1,p3.p5,1)>0 ) multiply 1/(1-1);
        endif;        
        .sort
        
************************************************************
*
* treat the scalar products
*
        if ( count(inte4,1) );
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        endif;        
        .sort

        if ( count(inte4,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        endif;        
        #call ACCU(E4_0)

        if ( count(inte4,1) );        
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        endif;        
        .sort

        #message numerator
        if ( count(inte4,1) );
        id  p1.p2 = 1/2 * (-1/s3m + 1/s1m + 1/s2m-M^2);
        id  p1.p3 = 1/2 * ( 1/s2m - 1/s1m - 1/s3m+M^2);
        id  p2.p3 = 1/2 * (-1/s1m + 1/s2m + 1/s3m-M^2);
        endif;        
        #call ACCU(E4_1)
        if ( count(inte4,1) );
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        endif;        
        .sort

        if ( count(inte4,1) );        
        id  p2.p4 = 1/2 * (-p5.p5 + 1/s2m + 1/s4m-2*M^2);
        id  p2.p5 = 1/2 * ( 1/s4m - 1/s2m - p5.p5);
        id  p4.p5 = 1/2 * (-1/s2m + 1/s4m + p5.p5);
        endif;        
        #call ACCU(E4_2)

        if ( count(inte4,1) );        
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        endif;        
        .sort

        if ( count(inte4,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        id p4.p4 = 1/s4m - M^2;
        endif;        
        .sort

* #include matad.info # time
* #include matad.info # print

************************************************************

* do recursion

        #message do recursion

        if ( count(inte4,1) );
        if (count(s4m,1)<=0) discard;
        if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
        if ( count(s1m,1) < count(s3m,1) ) multiply replace_(s1m,s3m,s3m,s1m);
        endif;        
        .sort

        #do i=1,1
                if ( count(inte4,1) );
                if ( match(s1m*s2m*s3m*s4m/p5.p5) > 0 );
                id s1m^n1? * s2m^n2? * s3m^n3? * s4m^n4? / p5.p5^n5? =
                s1m^n1 * s2m^n2 * s3m^n3 * s4m^n4 * (1/p5.p5)^n5 *
                deno(4-2*n5-n4,-2) * n4 * s4m * ( p5.p5 - 1/s2m )
                ;

                redefine i "0";
                
                endif;
                
                if (count(s4m,1)<=0) discard;
                if ( (count(s1m,1)<=0) && (count(s2m,1)<=0) ) discard;
                endif;        
                #call ACCU(E3)
                
        #enddo

* #include expandnomdeno

************************************************************

* identify simple integrals

        if ( count(inte4,1) );
        if ( (count(s1m,1) >0) && (count(s2m,1)<=0) && (count(s3m,1) >0) && 
        (count(s4m,1) >0) );
        id 1/s2m=p2.p2+M^2; 
        id p1=-p1;
        multiply replace_(p1,p3,p2,p1,p3,p6,s1m,s3m,s3m,s6m);
        multiply, intbn1/inte4;
        elseif ( (count(s1m,1) >0) && (count(s2m,1) >0) && (count(s3m,1) >0) && 
        (count(s4m,1) >0) && (count(p5.p5,1)>=0) );
        id p5=p4-p2;
        id p2=-p2;
        id p3=-p3;
        id p4=-p4;
        multiply replace_(p1,p5,p3,p1,s1m,s5m,s3m,s1m);
        multiply, intm5/inte4;
        elseif ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2;
        id p1=p2-p3;
        id p5=-p5;
        multiply replace_(p2,p4,p3,p5,p4,p6,p5,p1,s2m,s4m,s3m,s5m,s4m,s6m);
        multiply, intbm/inte4;
        elseif ( (count(s3m,1)<=0) );
        id 1/s3m=p3.p3+M^2;
        id p3=p2-p1;
        id p1=-p1;
        id p5=-p5;
        multiply replace_(p1,p5,p2,p4,p4,p6,p5,p1,s1m,s5m,s2m,s4m,s4m,s6m);
        multiply, intbm/inte4;
        else;
        multiply 1/(1-1);
        endif;
        endif;
        #call ACCU(topE4)

        #message - done
        
#endprocedure




************************************************************

#procedure mltadE3(s1m,s2m,s3m,p4,p5)
        
* discard massless tadpoles
        if ( count(inte3,1) );
        if ( (count(`s1m',1)<=0) && (count(`s2m',1)<=0) ) discard;
        if ( (count(`s3m',1)<=0) && (count(`s2m',1)<=0) ) discard;
        if ( (count(`p4'.`p4',1)>=0) ) discard;
        if ( (count(`p5'.`p5',1)>=0) ) discard;
        endif;        
        .sort

#endprocedure

#procedure redE3n2(s1m,s2m,s3m,p4,p5)

* reduce n2>0 to n2=1

        if ( count(inte3,1) );
        if (match(`s1m'*`s2m'^2*`s3m'/`p4'.`p4'/`p5'.`p5')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5 * 
        (-1)/M^2/(n2-1)*(
        + n1*M^2*`s1m'*`s2m'^-1
        - n1*`s2m'^-1
        + n3*M^2*`s2m'^-1*`s3m'
        - n3*`s2m'^-1
        - n4*`s2m'^-1
        - n5*`s2m'^-1
        - (n2-1)*`s2m'^-1
        + 3/2*`s2m'^-1*n
        );

        redefine i "0";

        endif;
        endif;
#endprocedure

#procedure redE3n45(s1m,s2m,s3m,p4,p5)

* reduce massless lines 4 and 5 to n4=1, n5=1;

* sort: n4 > n5
        if ( count(inte3,1) );
        if ( count(`p4'.`p4',1) > count(`p5'.`p5',1) )
        multiply replace_(`p4',`p5',`p5',`p4');
        endif;
        .sort

        if ( count(inte3,1) );        
        if (match(`s1m'*`s2m'*`s3m'/`p4'.`p4'^2/`p5'.`p5')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5 * 
        (-1)/M^2/(n4-1)*(
        - (n4-1)*`s2m'^-1
        + 2*`p4'.`p4'*n5
        + `p4'.`p4'*(n4-1)
        - `p4'.`p4'*n
        + `p5'.`p5'*(n4-1)
        );

        redefine i "0";

        endif;
        endif;
#endprocedure

#procedure redE3n13(s1m,s2m,s3m,p4,p5)

* reduce lines 1 and 3 to n1=1, n3=1;

* sort: n3 > n1
        if ( count(inte3,1) );
        if ( count(`s1m',1) > count(`s3m',1) )
        multiply replace_(`s1m',`s3m',`s3m',`s1m');
        endif;
        .sort

        if ( count(inte3,1) );        
        if (match(`s1m'*`s2m'*`s3m'^2/`p4'.`p4'/`p5'.`p5')>0);

        id `s1m'^n1?*`s2m'^n2?*`s3m'^n3?/`p4'.`p4'^n4?/`p5'.`p5'^n5? =
        `s1m'^n1*`s2m'^n2*`s3m'^n3/`p4'.`p4'^n4/`p5'.`p5'^n5 * 
        (-1)/M^2/(n3-1)*(
        + 2/3*n1*`s1m'*`s2m'^-1*`s3m'^-1
        - 2/3*n1*`s1m'*`s3m'^-2
        + 1/3*(n3-1)*`s1m'^-1
        - 1/3*(n3-1)*`s2m'^-1
        - (n3-1)*`s3m'^-1
        + 1/3*`s3m'^-1*n
        );

        redefine i "0";

        endif;
        endif;
#endprocedure



#procedure tope3
*
* this is topE3
*
        #message this is topE3

* As this topology results form other ones by shrinking one line
* in principle no scalar products may appear. However, one has to
* take care, that, e.g., in the topology DM the index of line 6 is 
* really reduced to 0 and not to negative values. Otherwise unwanted 
* scalar products appear.


************************************************************

* treat the scalar products

        #message numerator

        if ( count(inte3,1) );
        if ( count(p1.p5,1,p1.p4,1,p3.p4,1,p3.p5,1)>0 ) multiply 1/(1-1);
        endif;
        .sort

        if ( count(inte3,1) );
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        endif;        
        #call ACCU(E3_0)

        #call mltadE3(s1m,s2m,s3m,p5,p5)

        if ( count(inte3,1) );        
        id  p1.p2 = 1/2 * (-1/s3m + 1/s1m + p2.p2);
        id  p1.p3 = 1/2 * ( 1/s2m - 1/s1m - p3.p3);
        id  p2.p3 = 1/2 * (-1/s1m + 1/s2m + p3.p3);
        endif;        
        #call ACCU(E3_1)

        if ( count(inte3,1) );        
        id  p2.p4 = 1/2 * (-p5.p5 + p2.p2 + p4.p4);
        id  p2.p5 = 1/2 * ( p4.p4 - p2.p2 - p5.p5);
        id  p4.p5 = 1/2 * (-p2.p2 + p4.p4 + p5.p5);
        endif;        
        #call ACCU(E3_2)

        if ( count(inte3,1) );        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p3.p3 = 1/s3m - M^2;
        endif;        
        #call ACCU(E3_3)

        #call mltadE3(s1m,s2m,s3m,p5,p5)

* #include matad.info # time
* #include matad.info # print

************************************************************

* do recursion

        #message do recursion

* use symmetry: n1 < n3
*               n5 < n4

        if ( count(inte3,1) );        
        if ( count(s1m,1) > count(s3m,1) ) multiply replace_(s1m,s3m,s3m,s1m);
        if ( count(p4.p4,1) > count(p5.p5,1) ) multiply replace_(p4,p5,p5,p4);
        endif;        
        .sort

        #do i=1,1
                #call redE3n2(s1m,s2m,s3m,p4,p5)
                id n = num(4-2*ep);
                #call ACCU(E3)
                .sort
        #enddo

        #call mltadE3(s1m,s2m,s3m,p5,p5)

        #do i=1,1
                #call redE3n45(s1m,s2m,s3m,p4,p5)
                id n = num(4-2*ep);
                #call ACCU(E3)
                .sort
        #enddo

        #call mltadE3(s1m,s2m,s3m,p5,p5)

        id n = num(4-2*ep);
        #call ACCU(E3)

        #do i=1,1
                #call redE3n13(s1m,s2m,s3m,p4,p5)
                id n = num(4-2*ep);
                #call ACCU(E3)
                .sort
        #enddo

        #call mltadE3(s1m,s2m,s3m,p5,p5)

        id n = num(4-2*ep);
        #call ACCU(E3)
        .sort

************************************************************

* identify simple integrals

        if ( count(inte3,1) );
        if ( (count(s2m,1)<=0) );
        id 1/s2m=p2.p2+M^2; 
        id p2=p1+p3;
        id p1=-p1;
        multiply replace_(p1,p4,p3,p6,p4,p3,s1m,s4m,s3m,s6m);
        multiply, intbn2/inte3;

        elseif ( (count(s1m,1)<=0) );
        id 1/s1m=p1.p1+M^2;
        id p1=p2-p3;
        id p2=-p2;
        id p3=-p3;
        id p4=-p4;
        id p5=-p5;
        multiply replace_(p2,p6,p3,p5,p5,p1,s2m,s6m,s3m,s5m);
        multiply, intbm1/inte3;

        elseif ( (count(s3m,1)<=0) );
        id 1/s3m=p3.p3+M^2;
        id p3=p2-p1;
        id p2=-p2;
        id p4=-p4;
        id p5=-p5;
        multiply replace_(p1,p5,p2,p6,p5,p1,s1m,s5m,s2m,s6m);
        multiply, intbm1/inte3;

        elseif ( (count(s1m,1)==1)   && (count(s2m,1)==1)   && (count(s3m,1)==1) && 
        (count(p4.p4,1)==-1) && (count(p5.p5,1)==-1) );
        id s1m*s2m*s3m/p4.p4/p5.p5 = M^2 * miE3;
        Multiply int0/inte3;

        else;
        multiply 1/(1-1);
        endif;
        endif;

        #call ACCU(E3)

        #message - done
        
#endprocedure        



*************************************************
*
*      Recursion for simpler integrals
* 
*************************************************

#procedure symBNnom(p1,p2,p3,p4,p5,p6,x3,x4,x5,x6)
*
* sort: n6>=n3,n4,n5
*
        if ( count(intbn,1) );                
        if ( (count(`x3',1) > count(`x6',1)) && (count(`x3',1) > count(`x4',1))
        && (count(`x3',1) > count(`x5',1)) );
        id `p2'=-`p2';
        multiply replace_(`x3',`x6',`x6',`x4',`x4',`x5',`x5',`x3',
        `p3',`p6',`p6',`p4',`p4',`p5',`p5',`p3',
        `p2',`p1',`p1',`p2');
        endif;
        if ( (count(`x4',1) > count(`x6',1)) && (count(`x4',1) > count(`x5',1))
        && (count(`x4',1) > count(`x3',1)) );
        id `p1'=-`p1';
        multiply replace_(`x4',`x6',`x6',`x3',`x3',`x5',`x5',`x4',
        `p4',`p6',`p6',`p3',`p3',`p5',`p5',`p4',
        `p2',`p1',`p1',`p2');
        endif;
        if ( (count(`x5',1) > count(`x6',1)) && (count(`x5',1) > count(`x4',1))
        && (count(`x5',1) > count(`x3',1)) );
        id `p1'=-`p1';
        id `p2'=-`p2';
        multiply replace_(`x5',`x6',`x6',`x5',`x3',`x4',`x4',`x3',
        `p5',`p6',`p6',`p5',`p3',`p4',`p4',`p3');
        endif;
*
* sort: n4>n3 or, if n3=n4: n1>=n2
*
        if (count(`x3',1) > count(`x4',1));
        id `p1'=-`p1';
        id `p2'=-`p2';
        multiply replace_(`x3',`x4',`x4',`x3',
        `p3',`p4',`p4',`p3',`p2',`p1',`p1',`p2');
        endif;
        if ( (count(`x3',1) == count(`x4',1)) 
        && (count(`p1'.`p1',1) > count(`p2'.`p2',1)) );
        id `p1'=-`p1';
        id `p2'=-`p2';
        multiply replace_(`x3',`x4',`x4',`x3',
        `p3',`p4',`p4',`p3',`p2',`p1',`p1',`p2');
        endif;
        endif;
#endprocedure


#procedure symBN (p1,p2,x3,x4,x5,x6)
*
* sort: n6>=n3,n4,n5
*       n4>=n3
*
        if ( count(intbn,1) );        
        if (count(`x3',1) > count(`x6',1)) 
        multiply replace_(`x3',`x6',`x6',`x3',`x4',`x5',`x5',`x4');
        if (count(`x4',1) > count(`x6',1)) 
        multiply replace_(`x4',`x6',`x6',`x4',`x3',`x5',`x5',`x3');
        if (count(`x3',1) > count(`x4',1)) 
        multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
        if (count(`x5',1) > count(`x6',1)) 
        multiply replace_(`x5',`x6',`x6',`x5',`p1',`p2',`p2',`p1');
        endif;        
#endprocedure

#procedure BNtoBM (p1,p2,p3,p4,p5,p6,x3,x4,x5,x6)
        
* Changes the notation of the BM's which result when reducing BN.
*
* n6>=1 (after "symBN.prc" is used)
* There are several possibilities for the other n's:
*
* n3=0: 
*      n5!=0 (otherwise the result is ==0)
*      n4=0 or n4!=0 possible
* n4=0:
*      n5!=0 (otherwise the result is ==0)
*      n3=0  (because of "symBN.prc")
*      (-> see case n3=0)
* n5=0:
*      n3!=0 and n4!=0 (otherwise the result is ==0)
*
* So there are two different cases needed: n3=0 or n5=0.
        if ( count(intbn,1) );
        id `x3'^n3?neg_=(`p3'.`p3'+M^2)^(-n3);
        id `x4'^n4?neg_=(`p4'.`p4'+M^2)^(-n4);
        id `x5'^n5?neg_=(`p5'.`p5'+M^2)^(-n5);
        id `x6'^n6?neg_=(`p6'.`p6'+M^2)^(-n6);
        endif;        
        .sort
        
        if ( count(intbn,1) );        
        if((count(`x3',1) = 0) &&
        (count(`x4',1) >= 0)&&(count(`x5',1) >= 0)&&(count(`x6',1) >= 0) 
        ); 
        id `p3'=`p6'-`p1';
        multiply replace_(`p2',`p3',`x6',`x5',`x5',`x4',`x4',`x6',
        `p6',`p5',`p5',`p4',`p4',`p6');
        id `p1'=-`p1';
        elseif((count(`x5',1) = 0)&&
        (count(`x3',1) >= 0)&&(count(`x4',1) >= 0)&&(count(`x6',1) >= 0) 
        ); 
        id `p5'=`p2'+`p3';
        multiply replace_(`p1',`p3',`p2',`p1',`x3',`x5',
        `p3',`p5');
        endif;
        endif;
#endprocedure;



#procedure nomBN
*
* Decomposition of the numerator for type BN
*
        
        #call symBNnom(p1,p2,p3,p4,p5,p6,x3,x4,x5,x6)
        .sort
        if ( count(intbn,1) );
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;        
        .sort
        
        #do i=1,10
                
                #message pi.pi
                if ( count(intbn,1) );                
                id,once  p3.p3 = 1/x3 - M^2;
                id,once  p4.p4 = 1/x4 - M^2;
                id,once  p5.p5 = 1/x5 - M^2;
                id,once  p6.p6 = 1/x6 - M^2;
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;                
                
                #call ACCU(pi.pi)
                
        #enddo

        if ( count(intbn,1) );        
        id  p3.p3 = 1/x3 - M^2;
        id  p4.p4 = 1/x4 - M^2;
        id  p5.p5 = 1/x5 - M^2;
        id  p6.p6 = 1/x6 - M^2;
        
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        
        #call ACCU(pi.pi)

        #do i=1,10

                #message p1

                if ( count(intbn,1) );        
                id,once  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
                id,once  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
                id,once  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
                id,once  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;        
                
                #call ACCU(p1)
                
        #enddo

        if ( count(intbn,1) );
        id  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
        id  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
        id  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
        id  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );

        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;

        #call ACCU(p1)

        #do i=1,10

                #message p1,p2

                if ( count(intbn,1) );        
                id,once  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
                id,once  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
                id,once  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
                id,once  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;        
                
                #call ACCU(p1 p2)
                
        #enddo
        
        if ( count(intbn,1) );
        id  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
        id  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
        id  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
        id  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
        
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        
        #call ACCU(p1 p2)
        
        #do i=1,10
                
                #message p2,p3
                
                if ( count(intbn,1) );        
                id,once  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
                id,once  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
                id,once  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
                id,once  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;
                #call ACCU(p2 p3)
                
        #enddo
        
        if ( count(intbn,1) );
        id  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
        id  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
        id  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
        id  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
        
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        
        #call ACCU(p2 p3)
        
        #do i=1,10
                
                #message p4,p5,p6
                
                if ( count(intbn,1) );
                id,once  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
                id,once  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
                id,once  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
                
                if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
                if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
                endif;
                
                #call ACCU(p4 p5 p6)
                
        #enddo
        
        if ( count(intbn,1) );
        id  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
        id  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
        id  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
        
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        
        #call ACCU(p4 p5 p6)
        
        #call symBN(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );
        if ( (count(x3,1)>0)  && (count(x4,1)>0) 
        && (count(x5,1)<=0) && (count(x6,1)<=0)
        && ( (count(p1.p1,1)>=0)    
        || (count(p2.p2,1)>=0) )
        ) discard;
        endif;
        .sort
        
* #include expandnomdeno
        
        #message numerator decomposition done (BN)
        
* #include matad.info # time
* #include matad.info # print
        
#endprocedure






#procedure redBNn135 (p1,p2,x3,x4,x5,x6)
        if ( count(intbn,1) );
        repeat;
                if ( (count(`p1'.`p1',1) < 0)
                && (count(`x3',1) >= 1)  &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *deno(4-2*n1-n4-n6,-2)
                *(
                n4*`x4'*( `p1'.`p1' - 1/`x5' )
                +n6*`x6'*( `p1'.`p1' - 1/`x3' )
                )
                ;
                endif;
        endrepeat;
        endif;
        #call Conv2exact
#endprocedure

#procedure redBNn146 (p1,p2,x3,x4,x5,x6)
        
* sort: n6 < n4 < n3,n5
***if (count(`p1'.`p1',1) == 0)
***                 multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
        if ( count(intbn,1) );
        repeat;
                if ( (count(`p1'.`p1',1) < 0)
                && (count(`x3',1) > 0)   &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) > 0) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *deno(4-2*n1-n5-n3,-2)
                *(
                n5*`x5'*( `p1'.`p1' - 1/`x4' )
                +n3*`x3'*( `p1'.`p1' - 1/`x6' )
                )
                ;
                endif;
        endrepeat;
        endif;
#endprocedure

#procedure redBNn1p (p1,p2,x3,x4,x5,x6)

        #do i=1,1
                if ( count(intbn,1) );
                if ( (count(`p1'.`p1',1) > 0)
                && (count(`x3',1) >= 1)  &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *2*M^2*deno(3*4-2*(n1+n2+n3+n4+n5+n6),-6)
                *(
                -2*nom(4,-2)+4*(n1+1)+n3+n4+n5+n6
                -(n3*`x3'/`x6' + n4*`x4'/`x5' + n5*`x5'/`x4' + n6*`x6'/`x3')
                )/`p1'.`p1'
                ;
                
                redefine i "0";
                
                endif;
                endif;        
                .sort
                #call Conv2exact
        #enddo
        
#endprocedure

#procedure redBNn236 (p1,p2,x3,x4,x5,x6)
        if ( count(intbn,1) );
        repeat;
                
*
* sort: n2 < n1
*
***if (count(`p1'.`p1',1) > count(`p2'.`p2',1)) 
***                 multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                
                if ( (count(`p2'.`p2',1) < 0)
                && (count(`x3',1) > 0)   &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) > 0) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *deno(4-2*n2-n4-n5,-2)
                *(
                n4*`x4'*( `p2'.`p2' - 1/`x6' )
                +n5*`x5'*( `p2'.`p2' - 1/`x3' )
                )
                ;
                endif;
                
        endrepeat;
        endif;
        #call Conv2exact
#endprocedure

#procedure redBNn245 (p1,p2,x3,x4,x5,x6)
        if ( count(intbn,1) );
        repeat;
                if ( (count(`p2'.`p2',1) < 0)
                && (count(`x3',1) >= 1)  &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *deno(4-2*n2-n3-n6,-2)
                *(
                n6*`x6'*( `p2'.`p2' - 1/`x4' )
                +n3*`x3'*( `p2'.`p2' - 1/`x5' )
                )
                ;
                endif;
                
        endrepeat;
        endif;
        #call Conv2exact
#endprocedure

#procedure redBNn2p (p1,p2,x3,x4,x5,x6)
        
        #do i=1,1
                if ( count(intbn,1) );
                if ( (count(`p2'.`p2',1) > 0)
                && (count(`x3',1) >= 1)  &&  (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  

                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *2*M^2*deno(3*4-2*(n1+n2+n3+n4+n5+n6),-6)
                *(
                -2*nom(4,-2)+4*(n2+1)+n3+n4+n5+n6
                -(n4*`x4'/`x6' + n3*`x3'/`x5' + n5*`x5'/`x3' + n6*`x6'/`x4')
                )/`p2'.`p2'
                ;

                redefine i "0";

                endif;
                endif;        
                .sort
                #call Conv2exact
        #enddo

#endprocedure

#procedure redBNn342 (p1,p2,x3,x4,x5,x6)
        
        if ( count(intbn,1) );        
        repeat;
                
                if((count(`x3',1) >= 1)      &&  (count(`x4',1) >= 1) &&
                (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1)  );
                
* sort: n4 >= n3
                
                if (count(`x3',1) > count(`x4',1)) 
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                
                
* Now do the reduction via eq. (N1) resp. (5)
                
                if ( (count(`x4',1) > 1) && (count(`p2'.`p2',1) < 0) ); 
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^(n4-1) *(-1) * `x5'^n5  * `x6'^n6
                *1/2/(n4-1)/M^2
                *(
                +nom(4-2*(n4-1)-n1-n6,-2)
                +n1/`p1'.`p1'*(1/`x5'-1/`x4')
                +n6*`x6'     *(`p2'.`p2'-1/`x4'+2*M^2)
                )
                ;
                
                endif;
                endif;
                
        endrepeat;
        endif;
        #call Conv2exact
#endprocedure

#procedure redBNn3456 (p1,p2,x3,x4,x5,x6)

* do-enddo loop or not?
* example: reduction of BN(1,2,1,3,2,3)
*          with do-enddo loop: ~ 80s
*          without: ~ 375s

        #do i = 1,1
                
                if ( count(intbn,1) );
                if (   (count(`p1'.`p1',1) == 0) && (count(`p2'.`p2',1) == 0)
                && (count(`x3',1) >= 3)      && (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)      && (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *(
                ( -2*nom(4,-2) + 4*(n3-1) )/4/M^2/(n3-1) /x3
                -1/4/M^2/(n3-2)/(n3-1) * dala/x3^2
                )
                ;
                
                redefine i "0";
                
                endif;
                endif;        
                .sort
                #call Conv2exact
        #enddo
        
#endprocedure

#procedure redBNn3456exp (p1,p2,x3,x4,x5,x6)

* expand in ep
        
        #do i=1,1
                
                if ( count(intbn,1) );
                if (   (count(`p1'.`p1',1) == 0) && (count(`p2'.`p2',1) == 0)
                && (count(`x3',1) >= 3)      && (count(`x4',1) >= 1)  
                && (count(`x5',1) >= 1)      && (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^n5 * `x6'^n6
                *(
                num( -2*4 + 4*(n3-1)+ 4*ep )/4/M^2/(n3-1) /x3
                -1/4/M^2/(n3-2)/(n3-1) * dala/x3^2
                )
                ;
                
                redefine i "0";
                
                endif;
* topbn        
                endif;
                
*         repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
                #call ACCU(BNn3456)
                #call Conv2exact
        #enddo
        
#endprocedure

#procedure redBNn5 (p1,p2,x3,x4,x5,x6)

        #do i=1,1
                
* sort: n6>=n3,n4,n5>=1
*       n4>=n3
                if ( count(intbn,1) );     
                if((count(`x3',1) > 0)  &&  (count(`x4',1) > 0) &&
                (count(`x5',1) > 0)  &&  (count(`x6',1) > 0) );
                
                if (count(`x3',1) > count(`x6',1)) 
                multiply replace_(`x3',`x6',`x6',`x3',`x4',`x5',`x5',`x4');
                if (count(`x4',1) > count(`x6',1)) 
                multiply replace_(`x4',`x6',`x6',`x4',`x3',`x5',`x5',`x3');
                if (count(`x5',1) > count(`x6',1)) 
                multiply replace_(`x5',`x6',`x6',`x5',`p1',`p2',`p2',`p1');
                if (count(`x3',1) > count(`x4',1)) 
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                
* Now do the reduction via eq. (N15) resp. (6)
                
                if (count(`x5',1)>1);
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                *`x3'^n3 * `x4'^n4  * `x5'^(n5-1) *(-1) * `x6'^n6
                *(-1)/(n5-1)
                *(
                ( n3*`x3'+n4*`x4'+n6*`x6' )*(-1)
                -1/M^2*num(6-n1-n2-n3-n4-(n5-1)-n6-3*ep)
                )
                ;
                
                redefine i "0";
                endif;
                
                
                endif;
                endif;
                
                #call ACCU(BNn5)
                #call Conv2exact
        #enddo

#endprocedure

#procedure redBNn34 (p1,p2,x3,x4,x5,x6)
        
* at this stage: n5=1
*                n1,n2>=0
*                n3,n4,n6>0
        
        #do i=1,1
                if ( count(intbn,1) );
                if((count(`x3',1) >= 1)      &&  (count(`x4',1) >= 1) &&
                (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1)  );
                
* sort: n4 >= n3
                
                if (count(`x3',1) > count(`x4',1)) 
                multiply replace_(`x3',`x4',`x4',`x3',`p1',`p2',`p2',`p1');
                
* Now do the reduction via eq. (N1) resp. (5)
                
                if((count(`x4',1) > 1)); 
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^(n4-1) *(-1) * `x5'^n5  * `x6'^n6
                *1/2/(n4-1)/M^2
                *(
                +num(4-2*(n4-1)-n1-n6-2*ep)
                +n1/`p1'.`p1'*(1/`x5'-1/`x4')
                +n6*`x6'     *(`p2'.`p2'-1/`x4'+2*M^2)
                )
                ;
               
                redefine i "0"; 
                endif;
               
                
                endif;
                endif;
                #call ACCU(BNn34)
                #call Conv2exact
        #enddo
        
#endprocedure

#procedure redBNn12 (p1,p2,x3,x4,x5,x6)
        
* at this stage: n3=n4=n5=1
*                n1,n2 <=>0
*                n6>0
        
        #do i=1,1
                if ( count(intbn,1) );                
                if ( (count(`p1'.`p1',1) > 0)
                && (count(`x3',1) = 1)       &&  (count(`x4',1) = 1)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *M^2*den( 
                + 2*n1 + 2*n2 + 3*n6 - 7+ep*6
                )        
                *(-1)
                *(
                + p1.p1^-1*p2.p2^-1*x3^-1*x6^-1 * (  - n2*M^-2 )
                
                + p1.p1^-1*p2.p2^-1*x5^-1*x6^-1 * ( n2*M^-2 )
                
                - p1.p1^-1*x3^-1*x6 * ( 2*n6 )
                
                + p1.p1^-1*x3^-1 * (  - n6*M^-2 + M^-2 )
                
                - p1.p1^-1*x4^-1*x5 * ( 2 )
                
                - p1.p1^-1*x4*x5^-1 * ( 2 )
                
                + p1.p1^-1*x6^-1 * 1/M^2*num( - n2- n6 + 3 -ep*2)
                
                - p1.p1^-1 * num( 4 - 8*n1 - 4*n6 -ep*8 )
                
                
                ) ;
* N19b
                
                redefine i "0";
                
                endif;
                endif;                
                .sort
                #call Conv2exact
        #enddo
        
* #include expandnomdeno
        
        #do i=1,1
                
                if ( count(intbn,1) );        
                if ( (count(`p1'.`p1',1)<0 )
                && (count(`x3',1) = 1)       &&  (count(`x4',1) = 1)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *den(
                + 12  - 8*n1 - 4*n6 - ep*8
                )
                *(-1)
                *(
                + x3^-1*x6 * ( 2*n6 )
                
                - x3^-1 * (  - n6*M^-2 + M^-2 )
                
                + x4^-1*x5 * ( 2 )
                
                + x4*x5^-1 * ( 2 )
                
                - x6^-1 * 1/M^2*num( - n2- n6 + 3 - ep*2 )
                
                - p1.p1 * 1/M^2*num( + 2*n1 + 2*n2 + 3*n6 - 9+ep*6 )
                
                - p2.p2^-1*x3^-1*x6^-1 * (  - n2*M^-2 )
                
                - p2.p2^-1*x5^-1*x6^-1 * ( n2*M^-2 )
                
                
                );
* N18b
                
                redefine i "0";
                
                endif;
                endif;        
                .sort
                #call Conv2exact
        #enddo
        
* #include expandnomdeno
        
        #do i=1,1
                
                if ( count(intbn,1) );        
                if ( (count(`p2'.`p2',1)>0 )
                && (count(`x3',1) = 1)       &&  (count(`x4',1) = 1)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1) );  
                
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *M^2*den(
                + 2*n1 + 2*n2 + 3*n6 - 7 + ep*6
                )
                *(-1)
                *(
                + p1.p1^-1*p2.p2^-1*x4^-1*x6^-1 * (  - n1*M^-2 )
                
                + p1.p1^-1*p2.p2^-1*x5^-1*x6^-1 * ( n1*M^-2 )
                
                - p2.p2^-1*x3^-1*x5 * ( 2 )
                
                - p2.p2^-1*x3*x5^-1 * ( 2 )
                
                - p2.p2^-1*x4^-1*x6 * ( 2*n6 )
                
                + p2.p2^-1*x4^-1 * (  - n6*M^-2 + M^-2 )
                
                + p2.p2^-1*x6^-1 * 1/M^2*num(  - n1 - n6 + 3 - ep*2 )
                
                - p2.p2^-1 * num( 4 - 8*n2 - 4*n6 - ep*8 )
                
                
                ) ;
* N19a
                
                redefine i "0";
                
                endif;
                endif;        
                .sort
                #call Conv2exact
        #enddo
        
* #include expandnomdeno
        
        #do i=1,1
                if ( count(intbn,1) );
                if ( (count(`p2'.`p2',1)<0 )
                && (count(`x3',1) = 1)       &&  (count(`x4',1) = 1)  
                && (count(`x5',1) = 1)       &&  (count(`x6',1) >= 1) );  
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *den(
                + 12 - 8*n2 - 4*n6 - ep*8
                )
                *(-1)
                *(
                + x3^-1*x5 * ( 2 )
                
                + x3*x5^-1 * ( 2 )
                
                + x4^-1*x6 * ( 2*n6 )
                
                - x4^-1 * (  - n6*M^-2 + M^-2 )
                
                - x6^-1 * 1/M^2*num(   - n1 - n6 + 3 - ep*2 )
                
                - p1.p1^-1*x4^-1*x6^-1 * (  - n1*M^-2 )
                
                - p1.p1^-1*x5^-1*x6^-1 * ( n1*M^-2 )
                
                - p2.p2 * 1/M^2*num(  + 2*n1 + 2*n2 + 3*n6 - 9 + ep*6 )
                
                
                );
* N18a
                
                redefine i "0";
                
                endif;
                endif;        
                .sort
                #call Conv2exact
        #enddo
        
* #include expandnomdeno
        
#endprocedure

#procedure redBNn6 (p1,p2,x3,x4,x5,x6)
        
* at this stage: n1=n2=0
*                n3=n4=n5=1
* 		 n6>0    (?)             
        
        
        #do i=1,1
                
                if ( count(intbn,1) );        
                #ifdef `TABINT'
                        if((count(`x3',1) = 1)      &&  (count(`x4',1) = 1) &&
                        (count(`x5',1) = 1)      &&  (count(`x6',1) > 10) &&
                        (count(`p1'.`p1',1) = 0) &&  (count(`p2'.`p2',1) = 0)  );
                #endif
                
                #ifndef `TABINT'
                        
* old version without TabBN:
                        
                        if((count(`x3',1) = 1)      &&  (count(`x4',1) = 1) &&
                        (count(`x5',1) = 1)      &&  (count(`x6',1) > 1) &&
                        (count(`p1'.`p1',1) = 0) &&  (count(`p2'.`p2',1) = 0)  );
                #endif
                
*
* do the reduction via eq. (N20)
*
                
                id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? 
                * `x3'^n3? * `x4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 
                * `x3'^n3 * `x4'^n4  * `x5'^n5  * `x6'^n6
                *1/8/M^2/(n6-1)*den(4-n6-2*ep)*den(n6-3+2*ep)*den(n6-3+2*ep)
                *(-1)
                *(
                3*M**2*`x3'*`x4'*num(3 - n6-2*ep)*num(-4 + n6+2*ep)/(`x5'*`x6') 
                +
                3*M**2*`x4'*`x5'*num(3 - n6-2*ep)*num(-4 + n6+2*ep)/(`x3'*`x6') 
                +
                M**2*(-1 + n6)*`x3'*num(-4 + n6+2*ep)*num(-3 + n6+2*ep)/`x5' 
                +
                M**2*(-1 + n6)*`x5'*num(-4 + n6+2*ep)*num(-3 + n6+2*ep)/`x3' 
                -
                (-2 + n6)*`x3'*num(-4 + n6+2*ep)*num(-3 + n6+2*ep)/(`x5'*`x6') 
                -
                (-2 + n6)*`x5'*num(-4 + n6+2*ep)*num(-3 + n6+2*ep)/(`x3'*`x6') 
                +
                M**2*(-1 + n6)*n6*`x6'*num(-4 + n6+2*ep)*num(-3 + n6+2*ep)/`x4' 
                +
                (-2 + n6)*num(-5 + n6+3*ep)*num(-4 + n6+2*ep)*num(-3 + n6+2*ep)/
                (M**2*`x4'*`x6') 
                +       
                num(-5 + n6+3*ep)*num(-4 + n6+2*ep)**2*num(-3 + n6+2*ep)/
                (M^2*`x6'^2) 
                -    
                (-1 + n6)*num(-4 + n6+2*ep)*num(-3 + n6+2*ep)*
                num(-6 + 2*n6+3*ep)/`x4' 
                -
                num(-4 + n6+2*ep)*num(-3 + n6+2*ep)*
                nom(-58 + 41*n6 - 7*n6^2,54 - 20*n6,-12)/`x6'
                )
                
                ;
* N20a
                
                redefine i "0";
                
                endif;
                endif;
                .sort
                
                if ( count(intbn,1) );
                if ( (count(`x3',1)<=0)&&(count(`x5',1)<=0) ) discard;
                if ( (count(`x3',1)<=0)&&(count(`x6',1)<=0) ) discard;
                if ( (count(`x4',1)<=0)&&(count(`x5',1)<=0) ) discard;
                if ( (count(`x4',1)<=0)&&(count(`x6',1)<=0) ) discard;
                endif;

                #call Conv2exact                
                
                #call ACCU(n6)
                
        #enddo
        
#endprocedure


#procedure symmetryBN
*
* use symmetry for BN
*
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );            
*
* sort: n6 >= n5 >= n4 >= n3
*
        if (  ( count(x3,1) > count(x6,1) ) && ( count(x3,1) >= count(x5,1) )
        && ( count(x3,1) >= count(x4,1) ) ) 
        multiply replace_(x3,x6,x6,x3);
        if (  ( count(x4,1) > count(x6,1) ) && ( count(x4,1) >= count(x5,1) )
        && ( count(x4,1) >= count(x3,1) ) ) 
        multiply replace_(x4,x6,x6,x4);
        if (  ( count(x5,1) > count(x6,1) ) && ( count(x5,1) >= count(x4,1) )
        && ( count(x5,1) >= count(x3,1) ) ) 
        multiply replace_(x5,x6,x6,x5);
        if ( ( count(x3,1) > count(x5,1) ) && ( count(x3,1) >= count(x4,1) ) )
        multiply replace_(x3,x5,x5,x3);
        if ( ( count(x4,1) > count(x5,1) ) && ( count(x4,1) >= count(x3,1) ) )
        multiply replace_(x4,x5,x5,x4);
        if ( count(x3,1) > count(x4,1) ) 
        multiply replace_(x3,x4,x4,x3);
        endif;
        endif;        
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) );
        if ( count(x3,1) > count(x5,1) )  multiply replace_(x3,x5,x5,x3);
        if ( count(x4,1) > count(x6,1) )  multiply replace_(x4,x6,x6,x4);
        if ( count(x5,1) > count(x6,1) )  multiply replace_(x3,x4,x4,x3,x5,x6,x6,x5);
        endif;
        endif;        
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p2.p2,1)==0) );
        if ( count(x3,1) > count(x6,1) )  multiply replace_(x3,x6,x6,x3);
        if ( count(x4,1) > count(x5,1) )  multiply replace_(x4,x5,x5,x4);
        if ( count(x5,1) > count(x6,1) )  multiply replace_(x3,x4,x4,x3,x5,x6,x6,x5);
        endif;
        endif;        
        .sort
        
#endprocedure        



#procedure BNd0
id,only dala^( 0 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(4))*acc(
      275/12 + 2*ep^-3 + 23/3*ep^-2 + 3*ep^-1*z2 + 35/2*ep^-1 + 89/3*ep*z3 + 3/
      2*ep*z4 + 105/4*ep*z2 + 9/4*ep*z2^2 - 189/8*ep - 3*ep^2*z3*z2 + 525/2*
      ep^2*z3 - 649/4*ep^2*z4 - 6/5*ep^2*z5 + 275/8*ep^2*z2 + 69/8*ep^2*z2^2
       + 16*ep^2*B4 - 14917/48*ep^2 - 2*z3 + 23/2*z2
);
id,only dala^( 0 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(5))*acc(
      5/3 - ep^-3 - 7/3*ep^-2 - 3/2*ep^-1*z2 - 3*ep^-1 - 49/3*ep*z3 - 3/4*ep*
      z4 - 9/2*ep*z2 - 9/8*ep*z2^2 + 29*ep + 3/2*ep^2*z3*z2 - 109*ep^2*z3 + 
      329/4*ep^2*z4 + 3/5*ep^2*z5 + 5/2*ep^2*z2 - 21/8*ep^2*z2^2 - 8*ep^2*B4
       + 413/3*ep^2 + z3 - 7/2*z2
);
id,only dala^( 0 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(6))*acc(
      1/3 + 1/3*ep^-3 + 1/3*ep^-2 + 1/2*ep^-1*z2 + 1/3*ep^-1 - 8/3*ep*z3 + 43/
      4*ep*z4 + 1/2*ep*z2 + 3/8*ep*z2^2 - ep*B4 + 1/3*ep - 39/32*ep^2*z3*z2 + 
      5151/64*ep^2*z3 + 1/16*ep^2*z3^2 + 9/64*ep^2*z4*z2 - 8599/128*ep^2*z4 - 
      39/80*ep^2*z5 + 1/16*ep^2*z6 - 439/256*ep^2*z2 + 411/256*ep^2*z2^2 + 9/
      128*ep^2*z2^3 - 1/16*ep^2*ggam + 13/2*ep^2*B4 - 143503/1536*ep^2 - 8/3*
      z3 + 1/2*z2
);
id,only dala^( 0 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(7))*acc(
       - 63/8*ep*z4 + 3/4*ep*B4 + 69/128*ep^2*z3*z2 - 15965/256*ep^2*z3 - 3/64
      *ep^2*z3^2 - 27/256*ep^2*z4*z2 + 29925/512*ep^2*z4 + 69/320*ep^2*z5 - 3/
      64*ep^2*z6 + 1701/1024*ep^2*z2 - 945/1024*ep^2*z2^2 - 27/512*ep^2*z2^3
       + 3/64*ep^2*ggam - 45/8*ep^2*B4 + 144015/2048*ep^2 + 7/4*z3
);
id,only dala^( 0 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(8))*acc(
       - 3/8 + 63/32*ep*z3 - 63/32*ep*z4 + 3/16*ep*B4 + 69/512*ep^2*z3*z2 - 
      13949/1024*ep^2*z3 - 3/256*ep^2*z3^2 - 27/1024*ep^2*z4*z2 + 11781/2048*
      ep^2*z4 + 69/1280*ep^2*z5 - 3/256*ep^2*z6 - 603/4096*ep^2*z2 - 945/4096*
      ep^2*z2^2 - 27/2048*ep^2*z2^3 + 3/256*ep^2*ggam - 9/16*ep^2*B4 + 140943/
      8192*ep^2 + 7/16*z3
);
        
#endprocedure        

#procedure BNd
id,only dala^( 1 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(5))*acc(
      104/3 + 8/3*ep^-2 + 32/3*ep^-1 - 316/3*ep*z3 + 126*ep*z4 + 16*ep*z2 - 12
      *ep*B4 + 320/3*ep - 69/8*ep^2*z3*z2 + 31703/48*ep^2*z3 + 3/4*ep^2*z3^2
       + 27/16*ep^2*z4*z2 - 15077/32*ep^2*z4 - 69/20*ep^2*z5 + 3/4*ep^2*z6 + 
      1627/64*ep^2*z2 + 1137/64*ep^2*z2^2 + 27/32*ep^2*z2^3 - 3/4*ep^2*ggam + 
      46*ep^2*B4 - 308141/384*ep^2 - 28*z3 + 4*z2
);
id,only dala^( 1 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(6))*acc(
       - 2/3 - 2/3*ep^-2 - 2/3*ep^-1 + 16/3*ep*z3 - 63/2*ep*z4 - ep*z2 + 3*ep*
      B4 - 2/3*ep + 69/32*ep^2*z3*z2 - 46871/192*ep^2*z3 - 3/16*ep^2*z3^2 - 27/
      64*ep^2*z4*z2 + 27173/128*ep^2*z4 + 69/80*ep^2*z5 - 3/16*ep^2*z6 + 1445/
      256*ep^2*z2 - 1137/256*ep^2*z2^2 - 27/128*ep^2*z2^3 + 3/16*ep^2*ggam - 
      41/2*ep^2*B4 + 431021/1536*ep^2 + 7*z3 - z2
);
id,only dala^( 1 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(7))*acc(
       - 1/2 + 21/8*ep*z3 + 63/8*ep*z4 - 3/4*ep*B4 - 69/128*ep^2*z3*z2 + 16637/
      256*ep^2*z3 + 3/64*ep^2*z3^2 + 27/256*ep^2*z4*z2 - 35973/512*ep^2*z4 - 
      69/320*ep^2*z5 + 3/64*ep^2*z6 - 2469/1024*ep^2*z2 + 945/1024*ep^2*z2^2
       + 27/512*ep^2*z2^3 - 3/64*ep^2*ggam + 27/4*ep^2*B4 - 145039/2048*ep^2
       - 7/4*z3
);
id,only dala^( 1 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(8))*acc(
      3/4 - 7/2*ep*z3 + 63/16*ep*z4 - 3/8*ep*B4 - 3/8*ep - 69/256*ep^2*z3*z2
       + 14957/512*ep^2*z3 + 3/128*ep^2*z3^2 + 27/512*ep^2*z4*z2 - 13797/1024*
      ep^2*z4 - 69/640*ep^2*z5 + 3/128*ep^2*z6 + 603/2048*ep^2*z2 + 945/2048*
      ep^2*z2^2 + 27/1024*ep^2*z2^3 - 3/128*ep^2*ggam + 21/16*ep^2*B4 - 140943/
      4096*ep^2 - 7/8*z3
);
id,only dala^( 1 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(9))*acc(
      27/32 - 1183/256*ep*z3 + 441/128*ep*z4 - 21/64*ep*B4 + 153/128*ep - 483/
      2048*ep^2*z3*z2 + 75803/4096*ep^2*z3 + 21/1024*ep^2*z3^2 + 189/4096*ep^2
      *z4*z2 - 39123/8192*ep^2*z4 - 483/5120*ep^2*z5 + 21/1024*ep^2*z6 + 8829/
      16384*ep^2*z2 + 6615/16384*ep^2*z2^2 + 189/8192*ep^2*z2^3 - 21/1024*ep^2
      *ggam + 123/256*ep^2*B4 - 983913/32768*ep^2 - 49/64*z3
);
id,only dala^( 2 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(6))*acc(
       - 2/3 + 16/3*ep^-1 - 49/2*ep*z3 + 189/2*ep*z4 + 8*ep*z2 - 9*ep*B4 + 16/
      3*ep - 207/32*ep^2*z3*z2 + 141541/192*ep^2*z3 + 9/16*ep^2*z3^2 + 81/64*
      ep^2*z4*z2 - 75663/128*ep^2*z4 - 207/80*ep^2*z5 + 9/16*ep^2*z6 - 5359/
      256*ep^2*z2 + 2835/256*ep^2*z2^2 + 81/128*ep^2*z2^3 - 9/16*ep^2*ggam + 
      57*ep^2*B4 - 1297159/1536*ep^2 - 21*z3
);
id,only dala^( 2 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(7))*acc(
      4 - 63/4*ep*z3 - 1/2*ep - 147/8*ep^2*z3 + 567/8*ep^2*z4 + 6*ep^2*z2 - 27/
      4*ep^2*B4 + 4*ep^2
);
id,only dala^( 2 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(8))*acc(
       - 7/8 + 203/64*ep*z3 - 189/32*ep*z4 + 9/16*ep*B4 + 83/32*ep + 207/512*
      ep^2*z3*z2 - 54503/1024*ep^2*z3 - 9/256*ep^2*z3^2 - 81/1024*ep^2*z4*z2
       + 60543/2048*ep^2*z4 + 207/1280*ep^2*z5 - 9/256*ep^2*z6 - 273/4096*ep^2
      *z2 - 2835/4096*ep^2*z2^2 - 81/2048*ep^2*z2^3 + 9/256*ep^2*ggam - 183/64
      *ep^2*B4 + 423725/8192*ep^2 + 21/16*z3
);
id,only dala^( 2 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(9))*acc(
       - 81/32 + 3353/256*ep*z3 - 1323/128*ep*z4 + 63/64*ep*B4 - 351/128*ep + 
      1449/2048*ep^2*z3*z2 - 246337/4096*ep^2*z3 - 63/1024*ep^2*z3^2 - 567/
      4096*ep^2*z4*z2 + 145593/8192*ep^2*z4 + 1449/5120*ep^2*z5 - 63/1024*ep^2
      *z6 - 26487/16384*ep^2*z2 - 19845/16384*ep^2*z2^2 - 567/8192*ep^2*z2^3
       + 63/1024*ep^2*ggam - 453/256*ep^2*B4 + 2990907/32768*ep^2 + 147/64*z3
);
id,only dala^( 2 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(10))*acc(
       - 37/8 + 14007/512*ep*z3 - 567/32*ep*z4 + 27/16*ep*B4 - 8119/768*ep + 
      621/512*ep^2*z3*z2 - 150009/2048*ep^2*z3 - 27/256*ep^2*z3^2 - 243/1024*
      ep^2*z4*z2 + 17199/2048*ep^2*z4 + 621/1280*ep^2*z5 - 27/256*ep^2*z6 - 
      13107/4096*ep^2*z2 - 8505/4096*ep^2*z2^2 - 243/2048*ep^2*z2^3 + 27/256*
      ep^2*ggam - 477/512*ep^2*B4 + 10943303/73728*ep^2 + 63/16*z3
);
id,only dala^( 3 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(7))*acc(
       - 53/2 + 1617/16*ep*z3 - 567/8*ep*z4 + 27/4*ep*B4 + 9/8*ep + 621/128*
      ep^2*z3*z2 - 112437/256*ep^2*z3 - 27/64*ep^2*z3^2 - 243/256*ep^2*z4*z2
       + 36477/512*ep^2*z4 + 621/320*ep^2*z5 - 27/64*ep^2*z6 - 25395/1024*ep^2
      *z2 - 8505/1024*ep^2*z2^2 - 243/512*ep^2*z2^3 + 27/64*ep^2*ggam - 117/16
      *ep^2*B4 + 1246599/2048*ep^2 + 63/4*z3
);
id,only dala^( 3 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(8))*acc(
       - 53/8 + 2373/64*ep*z3 - 567/32*ep*z4 + 27/16*ep*B4 - 627/32*ep + 621/
      512*ep^2*z3*z2 - 34821/1024*ep^2*z3 - 27/256*ep^2*z3^2 - 243/1024*ep^2*
      z4*z2 - 72387/2048*ep^2*z4 + 621/1280*ep^2*z5 - 27/256*ep^2*z6 - 25395/
      4096*ep^2*z2 - 8505/4096*ep^2*z2^2 - 243/2048*ep^2*z2^3 + 27/256*ep^2*
      ggam + 207/64*ep^2*B4 + 1253511/8192*ep^2 + 63/16*z3
);
id,only dala^( 3 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(9))*acc(
      95/24 - 2037/128*ep*z3 + 567/32*ep*z4 - 27/16*ep*B4 - 1801/576*ep - 621/
      512*ep^2*z3*z2 + 146331/1024*ep^2*z3 + 27/256*ep^2*z3^2 + 243/1024*ep^2*
      z4*z2 - 122661/2048*ep^2*z4 - 621/1280*ep^2*z5 + 27/256*ep^2*z6 + 9011/
      4096*ep^2*z2 + 8505/4096*ep^2*z2^2 + 243/2048*ep^2*z2^3 - 27/256*ep^2*
      ggam + 747/128*ep^2*B4 - 36981277/221184*ep^2 - 63/16*z3
);
id,only dala^( 3 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(10))*acc(
      37/2 - 13503/128*ep*z3 + 567/8*ep*z4 - 27/4*ep*B4 + 7231/192*ep - 621/
      128*ep^2*z3*z2 + 10251/32*ep^2*z3 + 27/64*ep^2*z3^2 + 243/256*ep^2*z4*z2
       - 26271/512*ep^2*z4 - 621/320*ep^2*z5 + 27/64*ep^2*z6 + 13107/1024*ep^2
      *z2 + 8505/1024*ep^2*z2^2 + 243/512*ep^2*z2^3 - 27/64*ep^2*ggam + 693/
      128*ep^2*B4 - 11138159/18432*ep^2 - 63/4*z3
);
id,only dala^( 3 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(11))*acc(
      25555/512 - 2592219/8192*ep*z3 + 384993/2048*ep*z4 - 18333/1024*ep*B4 + 
      1778323/12288*ep - 421659/32768*ep^2*z3*z2 + 38887947/65536*ep^2*z3 + 
      18333/16384*ep^2*z3^2 + 164997/65536*ep^2*z4*z2 + 3768093/131072*ep^2*z4
       - 421659/81920*ep^2*z5 + 18333/16384*ep^2*z6 + 9231429/262144*ep^2*z2
       + 5774895/262144*ep^2*z2^2 + 164997/131072*ep^2*z2^3 - 18333/16384*ep^2
      *ggam - 10971/8192*ep^2*B4 - 7085415833/4718592*ep^2 - 42777/1024*z3
);
id,only dala^( 4 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(8))*acc(
      201/2 - 16611/32*ep*z3 + 2835/8*ep*z4 - 135/4*ep*B4 + 8267/48*ep - 3105/
      128*ep^2*z3*z2 + 432699/256*ep^2*z3 + 135/64*ep^2*z3^2 + 1215/256*ep^2*
      z4*z2 - 150633/512*ep^2*z4 - 621/64*ep^2*z5 + 135/64*ep^2*z6 + 77823/
      1024*ep^2*z2 + 42525/1024*ep^2*z2^2 + 1215/512*ep^2*z2^3 - 135/64*ep^2*
      ggam + 981/32*ep^2*B4 - 56655259/18432*ep^2 - 315/4*z3
);
id,only dala^( 4 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(9))*acc(
      201/4 - 20391/64*ep*z3 + 2835/16*ep*z4 - 135/8*ep*B4 + 15503/96*ep - 
      3105/256*ep^2*z3*z2 + 233367/512*ep^2*z3 + 135/128*ep^2*z3^2 + 1215/512*
      ep^2*z4*z2 + 121527/1024*ep^2*z4 - 621/128*ep^2*z5 + 135/128*ep^2*z6 + 
      77823/2048*ep^2*z2 + 42525/2048*ep^2*z2^2 + 1215/1024*ep^2*z2^3 - 135/
      128*ep^2*ggam - 639/64*ep^2*B4 - 51893467/36864*ep^2 - 315/8*z3
);
id,only dala^( 4 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(10))*acc(
       - 4111/128 + 288183/2048*ep*z3 - 65205/512*ep*z4 + 3105/256*ep*B4 - 
      24271/3072*ep + 71415/8192*ep^2*z3*z2 - 15029415/16384*ep^2*z3 - 3105/
      4096*ep^2*z3^2 - 27945/16384*ep^2*z4*z2 + 10223199/32768*ep^2*z4 + 14283/
      4096*ep^2*z5 - 3105/4096*ep^2*z6 - 1396713/65536*ep^2*z2 - 978075/65536*
      ep^2*z2^2 - 27945/32768*ep^2*z2^3 + 3105/4096*ep^2*ggam - 62793/2048*
      ep^2*B4 + 1440019661/1179648*ep^2 + 7245/256*z3
);
id,only dala^( 4 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(11))*acc(
       - 127775/512 + 12618879/8192*ep*z3 - 1924965/2048*ep*z4 + 91665/1024*ep
      *B4 - 8278295/12288*ep + 2108295/32768*ep^2*z3*z2 - 215177487/65536*ep^2
      *z3 - 91665/16384*ep^2*z3^2 - 824985/65536*ep^2*z4*z2 + 5799087/131072*
      ep^2*z4 + 421659/16384*ep^2*z5 - 91665/16384*ep^2*z6 - 46157145/262144*
      ep^2*z2 - 28874475/262144*ep^2*z2^2 - 824985/131072*ep^2*z2^3 + 91665/
      16384*ep^2*ggam - 91809/8192*ep^2*B4 + 36109955197/4718592*ep^2 + 213885/
      1024*z3
);
id,only dala^( 4 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(12))*acc(
       - 9304913/10240 + 199689105/32768*ep*z3 - 27910575/8192*ep*z4 + 1329075/
      4096*ep*B4 - 3762879793/1228800*ep + 30568725/131072*ep^2*z3*z2 - 
      2101636785/262144*ep^2*z3 - 1329075/65536*ep^2*z3^2 - 11961675/262144*
      ep^2*z4*z2 - 1120092435/524288*ep^2*z4 + 6113745/65536*ep^2*z5 - 1329075/
      65536*ep^2*z6 - 3378245559/5242880*ep^2*z2 - 418658625/1048576*ep^2*z2^2
       - 11961675/524288*ep^2*z2^3 + 1329075/65536*ep^2*ggam + 5836545/32768*
      ep^2*B4 + 61207571383387/2359296000*ep^2 + 3101175/4096*z3
);
id,only dala^( 5 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(9))*acc(
       - 31629/32 + 2983365/512*ep*z3 - 467775/128*ep*z4 + 22275/64*ep*B4 - 
      623279/256*ep + 512325/2048*ep^2*z3*z2 - 57051285/4096*ep^2*z3 - 22275/
      1024*ep^2*z3^2 - 200475/4096*ep^2*z4*z2 + 7390845/8192*ep^2*z4 + 102465/
      1024*ep^2*z5 - 22275/1024*ep^2*z6 - 11661147/16384*ep^2*z2 - 7016625/
      16384*ep^2*z2^2 - 200475/8192*ep^2*z2^3 + 22275/1024*ep^2*ggam - 57915/
      512*ep^2*B4 + 991203343/32768*ep^2 + 51975/64*z3
);
id,only dala^( 5 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(10))*acc(
       - 94887/128 + 10197495/2048*ep*z3 - 1403325/512*ep*z4 + 66825/256*ep*B4
       - 2628933/1024*ep + 1536975/8192*ep^2*z3*z2 - 99553095/16384*ep^2*z3 - 
      66825/4096*ep^2*z3^2 - 601425/16384*ep^2*z4*z2 - 67640265/32768*ep^2*z4
       + 307395/4096*ep^2*z5 - 66825/4096*ep^2*z6 - 34983441/65536*ep^2*z2 - 
      21049875/65536*ep^2*z2^2 - 601425/32768*ep^2*z2^3 + 66825/4096*ep^2*ggam
       + 360855/2048*ep^2*B4 + 2734270893/131072*ep^2 + 155925/256*z3
);
id,only dala^( 5 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(11))*acc(
      1157529/2560 - 17562825/8192*ep*z3 + 3529575/2048*ep*z4 - 168075/1024*ep
      *B4 + 41807523/102400*ep - 3865725/32768*ep^2*z3*z2 + 733970985/65536*
      ep^2*z3 + 168075/16384*ep^2*z3^2 + 1512675/65536*ep^2*z4*z2 - 412024725/
      131072*ep^2*z4 - 773145/16384*ep^2*z5 + 168075/16384*ep^2*z6 + 412489647/
      1310720*ep^2*z2 + 52943625/262144*ep^2*z2^2 + 1512675/131072*ep^2*z2^3
       - 168075/16384*ep^2*ggam + 2557575/8192*ep^2*B4 - 1076567281819/
      65536000*ep^2 - 392175/1024*z3
);
id,only dala^( 5 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(12))*acc(
      27914739/5120 - 586662615/16384*ep*z3 + 83731725/4096*ep*z4 - 3987225/
      2048*ep*B4 + 3576781533/204800*ep - 91706175/65536*ep^2*z3*z2 + 
      7103666775/131072*ep^2*z3 + 3987225/32768*ep^2*z3^2 + 35885025/131072*
      ep^2*z4*z2 + 2467138905/262144*ep^2*z4 - 18341235/32768*ep^2*z5 + 
      3987225/32768*ep^2*z6 + 10134736677/2621440*ep^2*z2 + 1255975875/524288*
      ep^2*z2^2 + 35885025/262144*ep^2*z2^3 - 3987225/32768*ep^2*ggam - 
      12193335/16384*ep^2*B4 - 20803897639049/131072000*ep^2 - 9303525/2048*z3
      
);
id,only dala^( 5 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(13))*acc(
      257571657/10240 - 11565975015/65536*ep*z3 + 771701175/8192*ep*z4 - 
      36747675/4096*ep*B4 + 77284295253/819200*ep - 845196525/131072*ep^2*z3*
      z2 + 82356662835/524288*ep^2*z3 + 36747675/65536*ep^2*z3^2 + 330729075/
      262144*ep^2*z4*z2 + 49817042415/524288*ep^2*z4 - 169039305/65536*ep^2*z5
       + 36747675/65536*ep^2*z6 + 93635373951/5242880*ep^2*z2 + 11575517625/
      1048576*ep^2*z2^2 + 330729075/524288*ep^2*z2^3 - 36747675/65536*ep^2*
      ggam - 547125435/65536*ep^2*B4 - 179459708314467/262144000*ep^2 - 
      85744575/4096*z3
);
id,only dala^( 6 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(10))*acc(
      11063547/640 - 225827595/2048*ep*z3 + 33041925/512*ep*z4 - 1573425/256*
      ep*B4 + 1328814969/25600*ep - 36188775/8192*ep^2*z3*z2 + 3142122795/
      16384*ep^2*z3 + 1573425/4096*ep^2*z3^2 + 14160825/16384*ep^2*z4*z2 + 
      564672465/32768*ep^2*z4 - 7237755/4096*ep^2*z5 + 1573425/4096*ep^2*z6 + 
      4036144221/327680*ep^2*z2 + 495628875/65536*ep^2*z2^2 + 14160825/32768*
      ep^2*z2^3 - 1573425/4096*ep^2*ggam - 2377755/2048*ep^2*B4 - 
      8361740207457/16384000*ep^2 - 3671325/256*z3
);
id,only dala^( 6 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(11))*acc(
      11063547/640 - 247855545/2048*ep*z3 + 33041925/512*ep*z4 - 1573425/256*
      ep*B4 + 1660721379/25600*ep - 36188775/8192*ep^2*z3*z2 + 1787157225/
      16384*ep^2*z3 + 1573425/4096*ep^2*z3^2 + 14160825/16384*ep^2*z4*z2 + 
      2150684865/32768*ep^2*z4 - 7237755/4096*ep^2*z5 + 1573425/4096*ep^2*z6
       + 4036144221/327680*ep^2*z2 + 495628875/65536*ep^2*z2^2 + 14160825/
      32768*ep^2*z2^3 - 1573425/4096*ep^2*ggam - 11818305/2048*ep^2*B4 - 
      7723909022337/16384000*ep^2 - 3671325/256*z3
);
id,only dala^( 6 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(12))*acc(
       - 25801737/2560 + 837975915/16384*ep*z3 - 77693175/2048*ep*z4 + 3699675/
      1024*ep*B4 - 2852820513/204800*ep + 85092525/32768*ep^2*z3*z2 - 
      29377113255/131072*ep^2*z3 - 3699675/16384*ep^2*z3^2 - 33297075/65536*
      ep^2*z4*z2 + 6737125185/131072*ep^2*z4 + 17018505/16384*ep^2*z5 - 
      3699675/16384*ep^2*z6 - 9327155391/1310720*ep^2*z2 - 1165397625/262144*
      ep^2*z2^2 - 33297075/131072*ep^2*z2^3 + 3699675/16384*ep^2*ggam - 
      84828465/16384*ep^2*B4 + 23395687784707/65536000*ep^2 + 8632575/1024*z3
);
id,only dala^( 6 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(13))*acc(
       - 1803001599/10240 + 79589911905/65536*ep*z3 - 5401908225/8192*ep*z4 + 
      257233725/4096*ep*B4 - 520384334211/819200*ep + 5916375675/131072*ep^2*
      z3*z2 - 669024439965/524288*ep^2*z3 - 257233725/65536*ep^2*z3^2 - 
      2315103525/262144*ep^2*z4*z2 - 299330421705/524288*ep^2*z4 + 1183275135/
      65536*ep^2*z5 - 257233725/65536*ep^2*z6 - 655447617657/5242880*ep^2*z2
       - 81028623375/1048576*ep^2*z2^2 - 2315103525/524288*ep^2*z2^3 + 
      257233725/65536*ep^2*ggam + 3241915245/65536*ep^2*B4 + 1280948932682229/
      262144000*ep^2 + 600212025/4096*z3
);
id,only dala^( 6 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(14))*acc(
       - 141289369929/143360 + 940453299765/131072*ep*z3 - 60455113425/16384*
      ep*z4 + 2878814925/8192*ep*B4 - 322042661882127/80281600*ep + 
      66212743275/262144*ep^2*z3*z2 - 4126734725715/1048576*ep^2*z3 - 
      2878814925/131072*ep^2*z3^2 - 25909334325/524288*ep^2*z4*z2 - 
      5140139914665/1048576*ep^2*z4 + 13242548655/131072*ep^2*z5 - 2878814925/
      131072*ep^2*z6 - 51380153918847/73400320*ep^2*z2 - 906826701375/2097152*
      ep^2*z2^2 - 25909334325/1048576*ep^2*z2^3 + 2878814925/131072*ep^2*ggam
       + 57593623185/131072*ep^2*B4 + 4608098483075964771/179830784000*ep^2 + 
      6717234825/8192*z3
);
id,only dala^( 7 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(11))*acc(
       - 298676151/640 + 12898079145/4096*ep*z3 - 893918025/512*ep*z4 + 
      42567525/256*ep*B4 - 82067986779/51200*ep + 979053075/8192*ep^2*z3*z2 - 
      127892119005/32768*ep^2*z3 - 42567525/4096*ep^2*z3^2 - 383107725/16384*
      ep^2*z4*z2 - 39719787345/32768*ep^2*z4 + 195810615/4096*ep^2*z5 - 
      42567525/4096*ep^2*z6 - 108704350593/327680*ep^2*z2 - 13408770375/65536*
      ep^2*z2^2 - 383107725/32768*ep^2*z2^3 + 42567525/4096*ep^2*ggam + 
      419645205/4096*ep^2*B4 + 216162350340381/16384000*ep^2 + 99324225/256*z3
      
);
id,only dala^( 7 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(12))*acc(
       - 298676151/512 + 69257958525/16384*ep*z3 - 4469590125/2048*ep*z4 + 
      212837625/1024*ep*B4 - 96404442027/40960*ep + 4895265375/32768*ep^2*z3*
      z2 - 329906695545/131072*ep^2*z3 - 212837625/16384*ep^2*z3^2 - 
      1915538625/65536*ep^2*z4*z2 - 370231197525/131072*ep^2*z4 + 979053075/
      16384*ep^2*z5 - 212837625/16384*ep^2*z6 - 108704350593/262144*ep^2*z2 - 
      67043851875/262144*ep^2*z2^2 - 1915538625/131072*ep^2*z2^3 + 212837625/
      16384*ep^2*ggam + 4141467225/16384*ep^2*B4 + 200405296878813/13107200*
      ep^2 + 496621125/1024*z3
);
id,only dala^( 7 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(13))*acc(
      11801595591/35840 - 57935155635/32768*ep*z3 + 5057200575/4096*ep*z4 - 
      240819075/2048*ep*B4 + 11646997128873/20070400*ep - 5538838725/65536*
      ep^2*z3*z2 + 1746535811265/262144*ep^2*z3 + 240819075/32768*ep^2*z3^2 + 
      2167371675/131072*ep^2*z4*z2 - 316504670265/262144*ep^2*z4 - 1107767745/
      32768*ep^2*z5 + 240819075/32768*ep^2*z6 + 4284570870513/18350080*ep^2*z2
       + 75858008625/524288*ep^2*z2^2 + 2167371675/262144*ep^2*z2^3 - 
      240819075/32768*ep^2*ggam + 4068936585/32768*ep^2*B4 - 
      514339430554699629/44957696000*ep^2 - 561911175/2048*z3
);
id,only dala^( 7 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(14))*acc(
      141289369929/17920 - 927018830115/16384*ep*z3 + 60455113425/2048*ep*z4
       - 2878814925/1024*ep*B4 + 312152405987097/10035200*ep - 66212743275/
      32768*ep^2*z3*z2 + 633398503185/16384*ep^2*z3 + 2878814925/16384*ep^2*
      z3^2 + 25909334325/65536*ep^2*z4*z2 + 4656499007265/131072*ep^2*z4 - 
      13242548655/16384*ep^2*z5 + 2878814925/16384*ep^2*z6 + 51380153918847/
      9175040*ep^2*z2 + 906826701375/262144*ep^2*z2^2 + 25909334325/131072*
      ep^2*z2^3 - 2878814925/16384*ep^2*ggam - 51835993335/16384*ep^2*B4 - 
      4698270428402960331/22478848000*ep^2 - 6717234825/1024*z3
);
id,only dala^( 7 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(15))*acc(
      238220921651049/4587520 - 1635454218850515/4194304*ep*z3 + 
      101922601597425/524288*ep*z4 - 4853457218925/262144*ep*B4 + 
      579727015161350697/2569011200*ep - 111629516035275/8388608*ep^2*z3*z2 + 
      1727778383527455/16777216*ep^2*z3 + 4853457218925/4194304*ep^2*z3^2 + 
      43681114970325/16777216*ep^2*z4*z2 + 10463116119841665/33554432*ep^2*z4
       - 22325903207055/4194304*ep^2*z5 + 4853457218925/4194304*ep^2*z6 + 
      86636809318439007/2348810240*ep^2*z2 + 1528839023961375/67108864*ep^2*
      z2^2 + 43681114970325/33554432*ep^2*z2^3 - 4853457218925/4194304*ep^2*
      ggam - 118494084664935/4194304*ep^2*B4 - 7432274218099299127131/
      5754585088000*ep^2 - 11324733510825/262144*z3
);
id,only dala^( 8 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(12))*acc(
      160848770193/8960 - 1036684721205/8192*ep*z3 + 68806683225/1024*ep*z4 - 3276508725/512*ep*B4 + 341786261184399/5017600*ep - 75359700675/16384*ep^2*z3*z2 + 6982233107535/65536*ep^2*z3 + 
      3276508725/8192*ep^2*z3^2 + 29488578525/32768*ep^2*z4*z2 + 4637475431505/65536*ep^2*z4 - 15071940135/8192*ep^2*z5 + 3276508725/8192*ep^2*z6 + 58509539860599/4587520*ep^2*z2 + 1032100248375/
      131072*ep^2*z2^2 + 29488578525/65536*ep^2*z2^3 - 3276508725/8192*ep^2*ggam - 51112404945/8192*ep^2*B4 - 5455732270292400027/11239424000*ep^2 - 7645187025/512*z3
);
id,only dala^( 8 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(13))*acc(
      482546310579/17920 - 3293538652215/16384*ep*z3 + 206420049675/2048*ep*z4 - 9829526175/1024*ep*B4 + 1160471750515317/10035200*ep - 226079102025/32768*ep^2*z3*z2 + 8506482668145/131072*ep^2*z3 + 
      9829526175/16384*ep^2*z3^2 + 88465735575/65536*ep^2*z4*z2 + 20517867884115/131072*ep^2*z4 - 45215820405/16384*ep^2*z5 + 9829526175/16384*ep^2*z6 + 175528619581797/9175040*ep^2*z2 + 
      3096300745125/262144*ep^2*z2^2 + 88465735575/131072*ep^2*z2^3 - 9829526175/16384*ep^2*ggam - 231973424235/16384*ep^2*B4 - 15218794973297619441/22478848000*ep^2 - 22935561075/1024*z3
);
id,only dala^( 8 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(14))*acc(
       - 17046569321181/1146880 + 87693448408335/1048576*ep*z3 - 7296490232325/131072*ep*z4 + 347451915825/65536*ep*B4 - 19853704100074653/642252800*ep + 7991394063975/2097152*ep^2*z3*z2 - 
      1153674051521355/4194304*ep^2*z3 - 347451915825/1048576*ep^2*z3^2 - 3127067242425/4194304*ep^2*z4*z2 + 308868717654315/8388608*ep^2*z4 + 1598278812795/1048576*ep^2*z5 - 347451915825/1048576*
      ep^2*z6 - 6196581969119883/587202560*ep^2*z2 - 109447353484875/16777216*ep^2*z2^2 - 3127067242425/8388608*ep^2*z2^3 + 347451915825/1048576*ep^2*ggam - 4111323438285/1048576*ep^2*B4 + 
      729927873089987876919/1438646272000*ep^2 + 810721136925/65536*z3
);
id,only dala^( 8 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(15))*acc(
       - 2143988294859441/4587520 + 14537892233481435/4194304*ep*z3 - 917303414376825/524288*ep*z4 + 43681114970325/262144*ep*B4 - 5084139420327568833/2569011200*ep + 1004665644317475/8388608*ep^2*z3
      *z2 - 22091822327149155/16777216*ep^2*z3 - 43681114970325/4194304*ep^2*z3^2 - 393130034732925/16777216*ep^2*z4*z2 - 87644998576339785/33554432*ep^2*z4 + 200933128863495/4194304*ep^2*z5 - 
      43681114970325/4194304*ep^2*z6 - 779731283865951063/2348810240*ep^2*z2 - 13759551215652375/67108864*ep^2*z2^2 - 393130034732925/33554432*ep^2*z2^3 + 43681114970325/4194304*ep^2*ggam + 
      988791446481615/4194304*ep^2*B4 + 68189056476855117705459/5754585088000*ep^2 + 101922601597425/262144*z3
);
id,only dala^( 8 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(16))*acc(
       - 9283383072342963/2621440
       + 458088928907232735/16777216*ep*z3
       - 27802673788949325/2097152*ep*z4
       + 1323936847092825/1048576*ep*B4
       - 23850251795116804299/1468006400*ep
       + 30450547483134975/33554432*ep^2*z3*z2
       - 38662985948860665/67108864*ep^2*z3
       - 1323936847092825/16777216*ep^2*z3^2
       - 11915431623835425/67108864*ep^2*z4*z2
       - 3284931390909449085/134217728*ep^2*z4
       + 6090109496626995/16777216*ep^2*z5
       - 1323936847092825/16777216*ep^2*z6
       - 3376277238051236709/1342177280*ep^2*z2
       - 417040106834239875/268435456*ep^2*z2^2
       - 11915431623835425/134217728*ep^2*z2^3
       + 1323936847092825/16777216*ep^2*ggam
       + 37451405023389315/16777216*ep^2*B4
       + 832015171392400061165131/9865003008000*ep^2
       + 3089185976549925/1048576*z3
      
);
id,only dala^( 9 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(13))*acc(
       - 53464091020587/57344
       + 1785557690498925/262144*ep*z3
       - 114365652951375/32768*ep*z4
       + 5445983473875/16384*ep*B4
       - 122808275214677739/32112640*ep
       + 125257619899125/524288*ep^2*z3*z2
       - 3680182482294465/1048576*ep^2*z3
       - 5445983473875/262144*ep^2*z3^2
       - 49013851264875/1048576*ep^2*z4*z2
       - 9956391706058175/2097152*ep^2*z4
       + 25051523979825/262144*ep^2*z5
       - 5445983473875/262144*ep^2*z6
       - 19445313496000941/29360128*ep^2*z2
       - 1715484794270625/4194304*ep^2*z2^2
       - 49013851264875/2097152*ep^2*z2^3
       + 5445983473875/262144*ep^2*ggam
       + 111720993348825/262144*ep^2*B4
       + 1735014661399909090017/71932313600*ep^2
       + 12707294772375/16384*z3
      
);
id,only dala^( 9 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(14))*acc(
       - 53464091020587/32768
       + 13108853982566475/1048576*ep*z3
       - 800559570659625/131072*ep*z4
       + 38121884317125/65536*ep*B4
       - 135639657059618619/18350080*ep
       + 876803339293875/2097152*ep^2*z3*z2
       - 4334585090074155/4194304*ep^2*z3
       - 38121884317125/1048576*ep^2*z3^2
       - 343096958854125/4194304*ep^2*z4*z2
       - 91652947309071225/8388608*ep^2*z4
       + 175360667858775/1048576*ep^2*z5
       - 38121884317125/1048576*ep^2*z6
       - 19445313496000941/16777216*ep^2*z2
       - 12008393559894375/16777216*ep^2*z2^2
       - 343096958854125/8388608*ep^2*z2^3
       + 38121884317125/1048576*ep^2*ggam
       + 1043454160187775/1048576*ep^2*B4
       + 1617118717193818460577/41104179200*ep^2
       + 88951063406625/65536*z3
      
);
id,only dala^( 9 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(15))*acc(
       + 813326473902609/917504
       - 21758397166032975/4194304*ep*z3
       + 1740083042872125/524288*ep*z4
       - 82861097279625/262144*ep*B4
       + 1071817104398323193/513802240*ep
       - 1905805237431375/8388608*ep^2*z3*z2
       + 252214205942836305/16777216*ep^2*z3
       + 82861097279625/4194304*ep^2*z3^2
       + 745749875516625/16777216*ep^2*z4*z2
       - 43237147387072275/33554432*ep^2*z4
       - 381161047486275/4194304*ep^2*z5
       + 82861097279625/4194304*ep^2*z6
       + 295759036854372087/469762048*ep^2*z2
       + 26101245643081875/67108864*ep^2*z2^2
       + 745749875516625/33554432*ep^2*z2^3
       - 82861097279625/4194304*ep^2*ggam
       + 618304316683725/4194304*ep^2*B4
       - 308062379950710456806731/10358253158400*ep^2
       - 193342560319125/262144*z3
      
);
id,only dala^( 9 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(16))*acc(
       + 9283383072342963/262144
       - 2265731156723764275/8388608*ep*z3
       + 139013368944746625/1048576*ep*z4
       - 6619684235464125/524288*ep*B4
       + 23330382343065598371/146800640*ep
       - 152252737415674875/16777216*ep^2*z3*z2
       + 1109492787558768795/33554432*ep^2*z3
       + 6619684235464125/8388608*ep^2*z3^2
       + 59577158119177125/33554432*ep^2*z4*z2
       + 15534971393300867025/67108864*ep^2*z4
       - 30450547483134975/8388608*ep^2*z5
       + 6619684235464125/8388608*ep^2*z6
       + 3376277238051236709/134217728*ep^2*z2
       + 2085200534171199375/134217728*ep^2*z2^2
       + 59577158119177125/67108864*ep^2*z2^3
       - 6619684235464125/8388608*ep^2*ggam
       - 176665530340203975/8388608*ep^2*B4
       - 848042540598718553654059/986500300800*ep^2
       - 15445929882749625/524288*z3
      
);
id,only dala^( 9 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(17))*acc(
       + 4455461991563872245/14680064
       - 160697427751201170375/67108864*ep*z3
       + 9531097044284300625/8388608*ep*z4
       - 453861764013538125/4194304*ep*B4
       + 2397090237098233593165/1644167168*ep
       - 10438820572311376875/134217728*ep^2*z3*z2
       - 15317182712726781075/33554432*ep^2*z3
       + 453861764013538125/67108864*ep^2*z3^2
       + 4084755876121843125/268435456*ep^2*z4*z2
       + 1257836303008199336625/536870912*ep^2*z4
       - 2087764114462275375/67108864*ep^2*z5
       + 453861764013538125/67108864*ep^2*z6
       + 1620417468151321066035/7516192768*ep^2*z2
       + 142966455664264509375/1073741824*ep^2*z2^2
       + 4084755876121843125/536870912*ep^2*z2^3
       - 453861764013538125/67108864*ep^2*ggam
       - 14406914497461640875/67108864*ep^2*B4
       - 76533513714741473633510089/11048803368960*ep^2
       - 1059010782698255625/4194304*z3
      
);
id,only dala^( 10 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(14))*acc(
       + 14415935810319315/229376
       - 496144255056274125/1048576*ep*z3
       + 30838155389724375/131072*ep*z4
       - 1468483589986875/65536*ep*B4
       + 7055073677791029711/25690112*ep
       - 33775122569698125/2097152*ep^2*z3*z2
       + 475866013268754675/4194304*ep^2*z3
       + 1468483589986875/1048576*ep^2*z3^2
       + 13216352309881875/4194304*ep^2*z4*z2
       + 3213069371906790375/8388608*ep^2*z4
       - 6755024513939625/1048576*ep^2*z5
       + 1468483589986875/1048576*ep^2*z6
       + 5243027333667327045/117440512*ep^2*z2
       + 462572330845865625/16777216*ep^2*z2^2
       + 13216352309881875/8388608*ep^2*z2^3
       - 1468483589986875/1048576*ep^2*ggam
       - 36415221368549625/1048576*ep^2*B4
       - 268744448432056156722479/172637552640*ep^2
       - 3426461709969375/65536*z3
      
);
id,only dala^( 10 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(15))*acc(
       + 14415935810319315/114688
       - 516703025316090375/524288*ep*z3
       + 30838155389724375/65536*ep*z4
       - 1468483589986875/32768*ep*B4
       + 7660542981824440941/12845056*ep
       - 33775122569698125/1048576*ep^2*z3*z2
       - 536700738631313025/4194304*ep^2*z3
       + 1468483589986875/524288*ep^2*z3^2
       + 13216352309881875/2097152*ep^2*z4*z2
       + 3953185101260175375/4194304*ep^2*z4
       - 6755024513939625/524288*ep^2*z5
       + 1468483589986875/524288*ep^2*z6
       + 5243027333667327045/58720256*ep^2*z2
       + 462572330845865625/8388608*ep^2*z2^2
       + 13216352309881875/4194304*ep^2*z2^3
       - 1468483589986875/524288*ep^2*ggam
       - 45226122908470875/524288*ep^2*B4
       - 250965662764022761850759/86318776320*ep^2
       - 3426461709969375/32768*z3
      
);
id,only dala^( 10 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(16))*acc(
       - 247744176316062345/3670016
       + 6853688262233323875/16777216*ep*z3
       - 529990823765143125/2097152*ep*z4
       + 25237658274530625/1048576*ep*B4
       - 71970861775037808129/411041792*ep
       + 580466140314204375/33554432*ep^2*z3*z2
       - 8804036862771489675/8388608*ep^2*z3
       - 25237658274530625/16777216*ep^2*z3^2
       - 227138924470775625/67108864*ep^2*z4*z2
       + 5012863848043324875/134217728*ep^2*z4
       + 116093228062840875/16777216*ep^2*z5
       - 25237658274530625/16777216*ep^2*z6
       - 90099261719123830335/1879048192*ep^2*z2
       - 7949862356477146875/268435456*ep^2*z2^2
       - 227138924470775625/134217728*ep^2*z2^3
       + 25237658274530625/16777216*ep^2*ggam
       - 91224023415107625/16777216*ep^2*B4
       + 6149083610771668125248173/2762200842240*ep^2
       + 58887869307238125/1048576*z3
      
);
id,only dala^( 10 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(17))*acc(
       - 49010081907202594695/14680064
       + 1750727532740040784125/67108864*ep*z3
       - 104842067487127306875/8388608*ep*z4
       + 4992479404148919375/4194304*ep*B4
       - 25868980865025415833375/1644167168*ep
       + 114827026295425145625/134217728*ep^2*z3*z2
       + 176280591928788013275/67108864*ep^2*z3
       - 4992479404148919375/67108864*ep^2*z3^2
       - 44932314637340274375/268435456*ep^2*z4*z2
       - 13226209122255997462875/536870912*ep^2*z4
       + 22965405259085029125/67108864*ep^2*z5
       - 4992479404148919375/67108864*ep^2*z6
       - 17824592149664531726385/7516192768*ep^2*z2
       - 1572631012306909603125/1073741824*ep^2*z2^2
       - 44932314637340274375/536870912*ep^2*z2^3
       + 4992479404148919375/67108864*ep^2*ggam
       + 151214271247861439625/67108864*ep^2*B4
       + 857977097255456339714679779/11048803368960*ep^2
       + 11649118609680811875/4194304*z3
      
);
id,only dala^( 10 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(18))*acc(
       - 736612970806862294715/23068672
       + 17254199182397282980875/67108864*ep*z3
       - 1002753839843001736875/8388608*ep*z4
       + 47750182849666749375/4194304*ep*B4
       - 4536240180767771393004135/28420603904*ep
       + 1098254205542335235625/134217728*ep^2*z3*z2
       + 52167136267654316872425/536870912*ep^2*z3
       - 47750182849666749375/67108864*ep^2*z3^2
       - 429751645647000744375/268435456*ep^2*z4*z2
       - 144843096640876362295875/536870912*ep^2*z4
       + 219650841108467047125/67108864*ep^2*z5
       - 47750182849666749375/67108864*ep^2*z6
       - 267900871146298726489245/11811160064*ep^2*z2
       - 15041307597645026053125/1073741824*ep^2*z2^2
       - 429751645647000744375/536870912*ep^2*z2^3
       + 47750182849666749375/67108864*ep^2*ggam
       + 1664634850495968495375/67108864*ep^2*B4
       + 1468117436873132388945941926583/2100851040583680*ep^2
       + 111417093315889081875/4194304*z3
      
);
id,only dala^( 11 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 1 )= M^(12-2*(15))*acc(
       - 4895022042320149755/917504
       + 172880232868742639625/4194304*ep*z3
       - 10471361223536049375/524288*ep*z4
       + 498636248739811875/262144*ep*B4
       - 2525482391978921265315/102760448*ep
       + 11468633721015673125/8388608*ep^2*z3*z2
       - 217000482835942575/2097152*ep^2*z3
       - 498636248739811875/4194304*ep^2*z3^2
       - 4487726238658306875/16777216*ep^2*z4*z2
       - 1249791802095111573375/33554432*ep^2*z4
       + 2293726744203134625/4194304*ep^2*z5
       - 498636248739811875/4194304*ep^2*z6
       - 1780289657253561679965/469762048*ep^2*z2
       - 157070418353040740625/67108864*ep^2*z2^2
       - 4487726238658306875/33554432*ep^2*z2^3
       + 498636248739811875/4194304*ep^2*ggam
       + 14255178523540849125/4194304*ep^2*B4
       + 29143571794243501409619277/230183403520*ep^2
       + 1163484580392894375/262144*z3
      
);
id,only dala^( 11 )*x3^( 1 )*x4^( 1 )*x5^( 1 )*x6^( 2 )= M^(12-2*(16))*acc(
       - 44055198380881347795/3670016
       + 1611769355677542686625/16777216*ep*z3
       - 94242251011824444375/2097152*ep*z4
       + 4487726238658306875/1048576*ep*B4
       - 24374068934029861705515/411041792*ep
       + 103217703489141058125/33554432*ep^2*z3*z2
       + 514734689915180952525/16777216*ep^2*z3
       - 4487726238658306875/16777216*ep^2*z3^2
       - 40389536147924761875/67108864*ep^2*z4*z2
       - 13258627573774925640375/134217728*ep^2*z4
       + 20643540697828211625/16777216*ep^2*z5
       - 4487726238658306875/16777216*ep^2*z6
       - 16022606915282055119685/1879048192*ep^2*z2
       - 1413633765177366665625/268435456*ep^2*z2^2
       - 40389536147924761875/134217728*ep^2*z2^3
       + 4487726238658306875/16777216*ep^2*ggam
       + 152231146651378612125/16777216*ep^2*B4
       + 245320904474093161783656693/920733614080*ep^2
       + 10471361223536049375/1048576*z3
      
);
id,only dala^( 11 )*x3^( 1 )*x4^( 1 )*x5^( 2 )*x6^( 2 )= M^(12-2*(17))*acc(
       + 128988185853912982515/20185088
       - 166983639811930470375/4194304*ep*z3
       + 25084817085899773125/1048576*ep*z4
       - 1194515099328560625/524288*ep*B4
       + 223161022164535726932765/12434014208*ep
       - 27473847284556894375/16777216*ep^2*z3*z2
       + 12218148059306990568075/134217728*ep^2*z3
       + 1194515099328560625/8388608*ep^2*z3^2
       + 10750635893957045625/33554432*ep^2*z4*z2
       + 107533950656601632625/67108864*ep^2*z4
       - 5494769456911378875/8388608*ep^2*z5
       + 1194515099328560625/8388608*ep^2*z6
       + 46911592014219542244645/10334765056*ep^2*z2
       + 376272256288496596875/134217728*ep^2*z2^2
       + 10750635893957045625/67108864*ep^2*z2^3
       - 1194515099328560625/8388608*ep^2*ggam
       + 106488897457721625/4194304*ep^2*B4
       - 126934479044733863730057397681/612748220170240*ep^2
       - 2787201898433308125/524288*z3
      
);
id,only dala^( 11 )*x3^( 1 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(18))*acc(
       + 2209838912420586884145/5767168
       - 51316929173928292615125/16777216*ep*z3
       + 3008261519529005210625/2097152*ep*z4
       - 143250548549000248125/1048576*ep*B4
       + 13381843747294800592240185/7105150976*ep
       - 3294762616627005706875/33554432*ep^2*z3*z2
       - 121993010438168384655525/134217728*ep^2*z3
       + 143250548549000248125/16777216*ep^2*z3^2
       + 1289254936941002233125/67108864*ep^2*z4*z2
       + 418485228485141059097625/134217728*ep^2*z4
       - 658952523325401141375/16777216*ep^2*z5
       + 143250548549000248125/16777216*ep^2*z6
       + 803702613438896179467735/2952790016*ep^2*z2
       + 45123922792935078159375/268435456*ep^2*z2^2
       + 1289254936941002233125/134217728*ep^2*z2^3
       - 143250548549000248125/16777216*ep^2*ggam
       - 4802903820089238488625/16777216*ep^2*B4
       - 1496060676386661860726847398183/175070920048640*ep^2
       - 334251279947667245625/1048576*z3
      
);
id,only dala^( 11 )*x3^( 2 )*x4^( 2 )*x5^( 2 )*x6^( 2 )= M^(12-2*(19))*acc(
       - 647054999806874812148743178689665/3670016
       + 26061617479066909827220400348552625/16777216*ep*z3
       - 1574474528171941139618272424998125/2097152*ep*z4
       + 74974977531997197124679639285625/1048576*ep*B4
       - 339300933656606691449481735675851835/411041792*ep
       + 1724424483235935533867631703569375/33554432*ep^2*z3*z2
       + 9426582702872406115821875150525025/134217728*ep^2*z3
       - 74974977531997197124679639285625/16777216*ep^2*z3^2
       - 674774797787974774122116753570625/67108864*ep^2*z4*z2
       - 190342828364736712461255010673785125/134217728*ep^2*z4
       + 344884896647187106773526340713875/16777216*ep^2*z5
       - 74974977531997197124679639285625/16777216*ep^2*z6
       - 199362554027182980342381272909017095/1879048192*ep^2*z2
       - 23617117922579117094274086374971875/268435456*ep^2*z2^2
       - 674774797787974774122116753570625/134217728*ep^2*z2^3
       + 74974977531997197124679639285625/16777216*ep^2*ggam
       + 2172267330046154842418614863676125/16777216*ep^2*B4
       + 4498668278196230427336746757612874762831/920733614080*ep^2
       + 174941614241326793290919158333125/1048576*z3
       - 108431217215972213061058560000*Mtep
      
);
        
#endprocedure        



#procedure reduceBNBN
*
* reduce BN to BM and simpler integrals
*

* Note: this file contains ALL procedures connected to the topology BN,
*       also the ones using the old reduction method.
 

* new method: 
* 1. Reduce indices of massless lines to zero:
*    n1 or n2 >0: "ordinary" triangle rule
*    n1 or n2 <0: "inverse" triangle rule
* 2. Use special recurrence relation for the resulting
*    BN(0,0,n3,n4,n5,n6) integrals.
*
*
* before each rec.rel.: use symmetry and expand (or use 'accun.prc')
*
************************************************************


************************************************************

* n1,n2 > 0

        #call symmetryBN
        if ( count(intbn,1) );
        if ( ( (count(x3,1)==count(x4,1)) || (count(x5,1)==count(x6,1)) )
        && (count(p1.p1,1) > count(p2.p2,1)) )
        multiply replace_(p1,p2,p2,p1);
        endif;        
        .sort
        
* #include expandnomdeno
        
* #include matad.info # time
* #include matad.info # print
        
* n1>0:
        
        if ( count(intbn,1) );        
        if ( count(x3,1,x5,1) < count(x4,1,x6,1) );
        #call redBNn135(p1,p2,x3,x4,x5,x6)
        else;
        #call redBNn146(p1,p2,x3,x4,x5,x6)
        endif;
        endif;        
        .sort
        
        #call symmetryBN
        
*      #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

* n2>0:
        if ( count(intbn,1) );
        if ( count(x4,1,x5,1) < count(x3,1,x6,1) );
        #call redBNn245(p1,p2,x3,x4,x5,x6)
        else;
        #call redBNn236(p1,p2,x3,x4,x5,x6)
        endif;
        endif;        
.sort

* n1,n2 < 0
        
        #call symmetryBN
        if ( count(intbn,1) );        
        if ( ( (count(x3,1)==count(x4,1)) || (count(x5,1)==count(x6,1)) )
        && (count(p1.p1,1) > count(p2.p2,1)) )
        multiply replace_(p1,p2,p2,p1);
        endif;        
        .sort

* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

* n1<0:

        #call redBNn1p(p1,p2,x3,x4,x5,x6)
        .sort
        
        #call symmetryBN
        
* #include expandnomdeno
        
************************************************************
        
* split epression

* * AFP
*         if ( match(x3*x4*x5*x6)>0 ) Multiply intd6^20;
*         b intd6,intbn;        
* * AFP        
*         Print+s;
*         .end        

        
* G diabnh7 = diabnbn;
* .sort
* skip diabnbn;
* if ( match(x3*x4*x5*x6)>0 ) discard;
* .sort
* g diabnbntmp1 = diabnbn-diabnh7;
* .sort
* drop diabnbn;
* .sort
* .store
* g diabnbn = diabnbntmp1;
* .sort

* #include matad.info # time
* #include matad.info # print

************************************************************

* n2<0:

        #call redBNn2p(p1,p2,x3,x4,x5,x6)
        .sort

* now: BN(0,0,n3,n4,n5,n6) or BM-type integrals

        #call symmetryBN
        
        if ( count(intbn,1) );        
        if ( ( (count(x3,1)==count(x4,1)) || (count(x5,1)==count(x6,1)) )
        && (count(p1.p1,1) > count(p2.p2,1)) )
        multiply replace_(p1,p2,p2,p1);
        endif;        
        .sort
        
* #include expandnomdeno

* If `EXPAND' is defined, `redBNn3456exp.prc' should be used,
* otherwise `redBNn3456.prc':
        
*         #call redBNn3456exp(p1,p2,x3,x4,x5,x6)
        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );
        multiply replace_(x3,x4,x4,x3);
        endif;
        endif;
*         #call redBNn3456exp(p1,p2,x3,x4,x5,x6)
        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );
        multiply replace_(x3,x5,x5,x3);
        endif;
        endif;
        
* #include expandnomdeno

*         #call redBNn3456exp(p1,p2,x3,x4,x5,x6)
        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );
        multiply replace_(x3,x6,x6,x3);
        endif;
        endif;
*         #call redBNn3456exp(p1,p2,x3,x4,x5,x6)
        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort

        #call symmetryBN
        
        if ( count(intbn,1) );        
        if ( ( (count(x3,1)==count(x4,1)) || (count(x5,1)==count(x6,1)) )
        && (count(p1.p1,1) > count(p2.p2,1)) )
        multiply replace_(p1,p2,p2,p1);
        endif;        
        .sort
        
* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

        .sort

************************************************************

* include table "BNdal.tbl"
        
        #message Use table for BN


*         Print+s;
*         .end
* Check if table is too small:
        
        if ( count(intbn,1) );       
*         if (count(dala,1) > 11) exit "table is too small"; 
        if (count(dala,1) > 11) multiply 1/(1-1); 

        #call BNd
        
* Replace temporarily x6 by s6m, in order not to touch
* terms, which are not listed in the previous table.

        if ( (count(dala,1)!=0) ) id x6=s6m;
        endif;        

        .sort

        if ( count(intbn,1) );                
        #call BNd0

        id s6m^n1? = x6^n1;
        endif;        
        .sort


*         if (count(x3,1,x4,1,x5,1,x6,1)>0) multiply 1/(1-1);

*         Print+s;
*         .end        

        
*         if ( count(intbn,1) );        
*         if (count(x3,1,x4,1,x5,1,x6,1)>0) multiply 1/(1-1);
*         endif;        

* `diabnbn' has to be zero, otherwise 
*  an entry is missing in the table

.sort

#endprocedure




*
*
* Old routine to generate tables and for cases
* which are not not tabulated 
*
* 
#procedure reduceBNnotab

        #message This is BN reduction without tables        
        
        #call symmetryBN
        
        #call ACCU{BN}
        
        #message 1
        
        #CALL redBNn5{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        
        #call ACCU{BN}
        

        #message 2
        
        #CALL redBNn34{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN
        
        #call ACCU{BN}
        
        #message 3
        .sort
        
        
*         #ifndef 'RECMAX'
*                 #define RECMAX "5"
*         #endif
*         #ifndef 'N6MAX'
*                 #define N6MAX "10"
*         #endif
*         #ifndef 'N6MIN'
*                 #define N6MIN "0"
*         #endif
*         #do i = 'N6MIN','N6MAX'
*                 #message n6 = 'i'
                
        
        #CALL redBNn12{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN


                
        
        #call ACCU{BN}
        
        if( count(intbn,1));
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;        
        .sort
        #message 5
        
        
        #CALL redBNn6{p1|p2|x3|x4|x5|x6}

        
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN

        
                
        #call ACCU{BN}
*                 if (count(acc,1)!=0) id accu(x?)=x;
        if( count(intbn,1));
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        .sort
        
                
        
        #call ACCU{BN}
*                 if (count(acc,1)!=0) id accun(x?)=x;
        #message 6
        
        
*         #enddo

*         Print+s;
*         .end        

        
        
        .sort
        
        #CALL redBNn12{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN
        
        
        #call ACCU{BN}
*         if (count(acc,1)!=0) id accun(x?)=x;
        
        #call ACCU(123)
        
        if( count(intbn,1));
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        .sort
        #message 7
        
        
        #CALL redBNn6{p1|p2|x3|x4|x5|x6}
        #call symBN{p1|p2|x3|x4|x5|x6}
        #call symmetryBN
        
        #message n6 done
        
        
        #call ACCU{BN}
        
        if( count(intbn,1));
        if ( (count(x3,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x3,1)<=0)&&(count(x6,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0)&&(count(x6,1)<=0) ) discard;
        endif;
        .sort
        #message 8
        
        #call symBN{p1,p2,x3,x4,x5,x6}
        
        #call BNtoBM{p1,p2,p3,p4,p5,p6,x3,x4,x5,x6}
        .sort
        
        #call ACCU(topBN)
        
        if ( count(intbn,1) );                
        if ( (count(x3,1)==0) ) multiply intbm/intbn;
        endif;

        if (count(intbm,1)==1);
        
        id  p4.p4 = 1/x4 - M^2;
        id  p5.p5 = 1/x5 - M^2;
        id  p6.p6 = 1/x6 - M^2;
        
        id  p1.p5 = 1/2 * ( p3.p3 + 1/x4  - p2.p2 - 1/x6);
        
        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
*
* sort n6>=n5>=n4
*  
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) ) 
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        
        if ( (count(x4,1)<=0) || (count(x5,1)<=0) || (count(x6,1)<=0) );
*
* sort n6>n5>n4
*
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) )
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;

        endif;
        .sort

*        id inttbl = 1;
        .sort

        #message - done
#endprocedure









#procedure topbn
*
* this is topbnbn
*
        #message this is topbnbn

        #message numerator
        #call nomBN

        #message do recursion

* Modified: applyed only to intbn        
*         #call reduceBNBN
*         id acc(x1?)*intbn = rat(x1,1)*int0;        
        #call reduceBNnotab
        

*         #call expansion(2)        
*         b ep;
*         Print+s;
*         .end        
        
        #call symBN{p1,p2,x3,x4,x5,x6}

        #call BNtoBM{p1,p2,p3,p4,p5,p6,x3,x4,x5,x6}
        .sort

        #call ACCU(topBN)

        if ( count(intbn,1) );                
        if ( (count(x3,1)==0)) multiply intbm/intbn;
        endif;
        
        if (count(intbm,1)==1);

        id  p4.p4 = 1/x4 - M^2;
        id  p5.p5 = 1/x5 - M^2;
        id  p6.p6 = 1/x6 - M^2;

        id  p1.p5 = 1/2 * ( p3.p3 + 1/x4  - p2.p2 - 1/x6);

        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
*
* sort n6>=n5>=n4
*  
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) ) 
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);

        if ( (count(x4,1)<=0) || (count(x5,1)<=0) || (count(x6,1)<=0) );
*
* sort n6>n5>n4
*
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) )
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;

        endif;
        .sort

        
*        id inttbl = 1;
        .sort

        #message - done

#endprocedure        



#procedure topbn2
*
* this is topbn2
*
        #-
        #message this is topbn2

        #message numerator

        if( count(intbn2,1));        
        id  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
        id  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
        endif;        
        #call ACCU(BN2 1)

        if( count(intbn2,1));        
        id  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
        id  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );
        endif;        
        #call ACCU(BN2 2)

        if( count(intbn2,1));        
        id  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
        id  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
        endif;        
        #call ACCU(BN2 3)

        if( count(intbn2,1));        
        id  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
        id  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
        endif;        
        #call ACCU(BN2 4)

        if( count(intbn2,1));        
        id  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
        id  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN2 5)

        if( count(intbn2,1));        
        id  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
        id  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN2 6)

        if( count(intbn2,1));        
        id  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
        id  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
        id  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
        endif;        
        #call ACCU(BN2 7)

*
* Warning! 
*
        if( count(intbn2,1));        
        id  1/x3 = M^2 + p3.p3;
        id  1/x5 = M^2 + p5.p5;
        endif;        
        #call ACCU(BN2 8)

        if( count(intbn2,1));        
        id  p6.p6 = 1/x6 - M^2;
        id  p4.p4 = 1/x4 - M^2;
        endif;        
        #call ACCU(BN2 7)

        #message do recursion

        if( count(intbn2,1));        
        id p5=-p5;
        multiply replace_(p6,p5,p4,p6,p1,p4,p3,p2,p5,p1,p2,p3,x6,x5,x4,x6);
        multiply intbm1/intbn2;
        endif;        
        .sort

        #message - done
        
#endprocedure        



************************************************************

#procedure symBN3 (p1,p2,p3,p4,p5,x6)
        if( count(intbn3,1));
        if (count(`p1'.`p1',1) >= count(`p2'.`p2',1))
        multiply replace_(`p1',`p2',`p2',`p1',`p4',`p3',`p3',`p4');
        endif;
#endprocedure



#procedure triown(P,PA,PB,P1,P3)
*
*   Routine solves the triangle recursion
*                     P
*      N1 -------------------------- N3
*          P1    \    N2    /    P3
*                 \        /
*                A \      / B
*                PA \    / PB
*                    \  /
*                     \/
*
*   A,B,N1,N2,N3 are the powers of the denominators
*       We assume that A and B are integers here. Otherwise see triangl2.
*   P is the momentum in N2. We need it here to determine the extra
*   momenta in the numerator.
*
*   We follow the algorithm of F.V.Tkachov Theor. Mat. Fiz. 56(1983)350.
*
*   There are three sums:
*   One in which the power of 1 becomes zero, one for 2 and one for 3.
*   Each sum has 4 constants to sum over, except that for each
*   one of them has its maximal value:
*   k1 = (a-A)-k2
*   k2 = N1-n1          sum1: N1
*   k3 = (b-B)-k4
*   k4 = N3-n3          sum3: N3
*   k1+k3 = N2-n2       sum2: N2
*   The gamma functions are evaluated using the tables in pochtabl.prc
*
*   It turns out that the easy cases are faster when the regular recursion
*   is used. We do that here is there are fewer than 6 powers in P,P1 and P3
*   combined. Test timings don't use these special cases.
*
*   Programmed by J.A.M.Vermaseren 28-oct-1990 + 13-nov-1990
*
*   Test versus benzbar:
*       1/p1.p1^3/p.p/p3.p3/pa.pa/pb.pb -> 0.96 sec versus 0.80 sec
*       1/p1.p1^3/p.p^3/p3.p3/pa.pa/pb.pb -> 3.07 sec versus 20.34 sec
*       1/p1.p1^3/p.p^3/p3.p3^3/pa.pa/pb.pb -> 6.21 sec versus 547.11 sec
*   triangle generates exactly the right number of terms
*   (times on Atari TT)
*
*   For the easy cases we take the original recursion, because this is much
*   faster. This is done for the combined powers of P1,P2,P3 at most 4
*   in the denominator. IT IS VERY DANGEROUS TO INCREASE THIS NUMBER!
*   In principle each two steps in the recursion can generate one
*   extra pole that has to be cancelled between the terms. This can
*   cause problems with the truncation of the powers of ep. We can
*   tolerate only one such pole here.
*
id  'P' = x*'P';
if ( count('P1'.'P1',1,'P3'.'P3',1,'P'.'P',1) < -4 );
id  x^m?/'P1'.'P1'^N1?/'P'.'P'^N2?/'P3'.'P3'^N3?/'PA'.'PA'^A?/'PB'.'PB'^B? =
        ftri(N1,N2,N3,A,B,m)/fac_(A-1)/fac_(B-1);
id  ftri(N1?,N2?,N3?,A?,B?,m?) =
        +sum_(k4,0,N3-1,(-1)^(N1+k4)*po(5+m-2*N2-A-B-N1-k4,-2)
            *sum_(k1,0,N2-1,sum_(k3,0,N2-k1-1,
            fac_(k1+N1+k3+k4-1)/fac_(k3)
            *poinv(5+m-2*N2-A-B+k1+k3,-2)
            *ftri(0,N2-k1-k3,N3-k4,A+k1+N1,B+k3+k4)
         )/fac_(k1))/fac_(k4))/fac_(N1-1)

        +sum_(k2,0,N1-1,sum_(k4,0,N3-1,(-1)^(k2+k4)*fac_(k2+N2+k4-1)
            *po(4+m-2*N2-A-B-k2-k4,-2)
            *sum_(k1,0,N2,1/fac_(k1)/fac_(N2-k1)
            *ftri(N1-k2,0,N3-k4,A+k1+k2,B+N2-k1+k4)
         )/fac_(k2)/fac_(k4)))*N2*poinv(4+m-2*N2-A-B+N2,-2)

        +sum_(k2,0,N1-1,(-1)^(k2+N3)*po(5+m-2*N2-A-B-k2-N3,-2)
            *sum_(k3,0,N2-1,sum_(k1,0,N2-k3-1,
            fac_(k1+k2+k3+N3-1)/fac_(k1)
            *poinv(5+m-2*N2-A-B+k1+k3,-2)
            *ftri(N1-k2,N2-k1-k3,0,A+k1+k2,B+k3+N3)
         )/fac_(k3))/fac_(k2))/fac_(N3-1)
    ;
id  ftri(N1?,N2?,N3?,A?,B?) = fac_(A-1)*fac_(B-1)
    /'P1'.'P1'^N1/'P'.'P'^N2/'P3'.'P3'^N3/'PA'.'PA'^A/'PB'.'PB'^B;
*
*   Evaluate the gamma functions.
*   We multiply numerator and denominator with gamma(1-2*ep)
*   Then po(x,-2) represents a normalized gamma function.
*   The following identities can then be applied:
*   (the table has overal factors ep taken out)
*
id  po(1,-2) = 1;
id  poinv(1,-2) = 1;

***id  po(x1?neg0_,-2) = -acc(PO(x1,-2))/2/ep;
***id  po(x1?,-2) = acc(PO(x1,-2));
***id  poinv(x1?neg0_,-2) = -2*acc(POINV(x1,-2))*ep;
***id  poinv(x1?,-2) = acc(POINV(x1,-2));

else;

  repeat;
    if ( match(1/'P'.'P'/'P1'.'P1'/'P3'.'P3') > 0 );
      id  x^x4?/'PA'.'PA'^x1?/'PB'.'PB'^x2?/'P'.'P'^x3? =
          x^x4 /'PA'.'PA'^x1 /'PB'.'PB'^x2 /'P'.'P'^x3 *(
            +x2*'P'.'P'/'PB'.'PB'-x2*'P3'.'P3'/'PB'.'PB'
            +x1*'P'.'P'/'PA'.'PA'-x1*'P1'.'P1'/'PA'.'PA'
          )*deno(4+x4-x1-x2-2*x3,-2);
***      id,many,deno(0,-2)   = -1/2/ep;
***      id,many,deno(x1?,-2) = accm(1-2*ep/x1)/x1;
    endif;
  endrepeat;
***  repeat id accm(x1?)*accm(x2?) = accm(x1*x2);
***  id  accm(x1?) = accm(x1-1);
***  id  accm(x1?) = accm(x1-x1^2,x1^3);
***  id  accm(x1?,x2?) = acc(1-x1-x2+x2*x1+x2^2);
  id  x = 1;

endif;

* repeat id acc(x1?)*acc(x2?) = acc(x1*x2);

#endprocedure




*--#[ ntriangle :
*
#procedure ntriangle(p,pa,pb,p1,p2)
*
*   Routine solves the triangle recursion
*                     p
*      n1 -------------------------- n2
*          p1    \    n0    /    p2
*                 \        /
*               n3 \      / n4
*                pa \    / pb
*                    \  /
*                     \/
*
*   n3,n4,n0,n1,n2 are the powers of the denominators
*   eppa and eppb are 1/pa.pa^ep and 1/pb.pb^ep
*   p is the momentum in n0. We need it here to determine the extra
*   momenta in the numerator.
*
id	`p' = xpower*`p';
if ( count(`p1'.`p1',1,`p2'.`p2',1,`p'.`p',1) >= -4 );
*
*	Here we can just use the recursion
*
	repeat;
		id	xpower^n8?/`p1'.`p1'^n1?pos_/`p2'.`p2'^n2?pos_/`p'.`p'^n?pos_/`pa'.`pa'^n3?/`pb'.`pb'^n4?
			*ep`pa'^x1?*ep`pb'^x2? = xpower^n8/`p1'.`p1'^n1/`p2'.`p2'^n2/`p'.`p'^n/`pa'.`pa'^n3/`pb'.`pb'^n4
				*ep`pa'^x1*ep`pb'^x2*(
					+num(n3+x1*ep)*(`p'.`p'/`pa'.`pa'-`p1'.`p1'/`pa'.`pa')
					+num(n4+x2*ep)*(`p'.`p'/`pb'.`pb'-`p2'.`p2'/`pb'.`pb')
				)*den(4-2*ep+n8-2*n-n3-n4-x1*ep-x2*ep);
	endrepeat;
	id	num(x?number_) = x;
	id	den(x?number_) = 1/x;
	id	num(x?)*den(x?) = 1;
	id	num(x?) = rat(x,1);
	id	den(x?) = rat(1,x);
else;
*
*	Here we have to use the general formula
*
	id	xpower^n8?/`p1'.`p1'^n1?pos_/`p2'.`p2'^n2?pos_/`p'.`p'^n?pos_/`pa'.`pa'^n3?/`pb'.`pb'^n4?
			*ep`pa'^x1?*ep`pb'^x2? = ftriangle(n8,n,n3,x1,n4,x2,n1,n2);
	id	ftriangle(n?,n0?,n3?,x1?,n4?,x2?,n1?,n2?) =
		+sum_(isum1,0,n2-1,sum_(isum2,0,n0-1,sum_(isum3,0,n0-1-isum2,
			ftriangle(n,n0-isum2-isum3,n3+n1+isum2,x1,n4+isum1+isum3,x2,0,n2-isum1)*sign_(n1+isum1)
				*fac_(n1+isum1+isum2+isum3-1)
				*Pochhammer(n1+isum2,n3+x1*ep)
				*Pochhammer(isum3+isum1,n4+x2*ep)
				*Pochhammer(-n1-isum1-isum2-isum3,isum2+isum3+5-2*ep+n-2*n0-n3-n4-x1*ep-x2*ep)
				*invfac_(isum1)*invfac_(isum2)*invfac_(isum3)*invfac_(n1-1)
		)))
		+sum_(isum1,0,n1-1,sum_(isum2,0,n2-1,sum_(isum3,0,n0,
			ftriangle(n,0,n3+isum1+isum3,x1,n4+n0+isum2-isum3,x2,n1-isum1,n2-isum2)*sign_(isum1+isum2)
				*fac_(n0+isum1+isum2-1)
				*Pochhammer(isum3+isum1,n3+x1*ep)
				*Pochhammer(n0-isum3+isum2,n4+x2*ep)
				*Pochhammer(-isum1-isum2-n0,4-2*ep+n-n0-n3-n4-x1*ep-x2*ep)
				*n0*invfac_(n0-isum3)
				*invfac_(isum1)*invfac_(isum2)*invfac_(isum3)
		)))
		+sum_(isum1,0,n1-1,sum_(isum2,0,n0-1,sum_(isum3,0,n0-1-isum2,
			ftriangle(n,n0-isum2-isum3,n3+isum1+isum3,x1,n4+n2+isum2,x2,n1-isum1,0)*sign_(n2+isum1)
				*fac_(n2+isum1+isum2+isum3-1)
				*Pochhammer(n2+isum2,n4+x2*ep)
				*Pochhammer(isum3+isum1,n3++x1*ep)
				*Pochhammer(-n2-isum1-isum2-isum3,isum2+isum3+5-2*ep+n-2*n0-n3-n4-x1*ep-x2*ep)
				*invfac_(isum1)*invfac_(isum2)*invfac_(isum3)*invfac_(n2-1)
		)));
	repeat id Pochhammer(n?pos_,x?) = Pochhammer(n-1,x)*num(n-1+x);
	repeat id Pochhammer(n?neg_,x?) = Pochhammer(n+1,x)*den(n+x);
	id	Pochhammer(0,x?) = 1;
	id	num(x?number_) = x;
	id	den(x?number_) = 1/x;
	id	num(x?)*den(x?) = 1;
	id	num(x?) = rat(x,1);
	id	den(x?) = rat(1,x);
	id	ftriangle(n?,n0?,n3?,x1?,n4?,x2?,n1?,n2?) =
		ep`pa'^x1*ep`pb'^x2/`p'.`p'^n0/`p1'.`p1'^n1/`p2'.`p2'^n2/`pa'.`pa'^n3/`pb'.`pb'^n4;
endif;
id	xpower = 1;
*
#endprocedure
*
*--#] ntriangle : 



#procedure topbn3
*
* this is topbn3
*
        #-
        
        #message this is topbn3

************************************************************

        #message numerator
        
        if( count(intbn3,1));
        id  p1.p2 = 1/2 * ( 1/x4 + 1/x3 - 1/x5 - 1/x6 );
        id  p1.p3 = 1/2 * ( 1/x6 - 1/x3 - p1.p1 );
        endif;
        #call ACCU(BN3 1)

        if( count(intbn3,1));
        id  p1.p4 = 1/2 * (-1/x5 + 1/x4 + p1.p1 );
        id  p1.p5 = 1/2 * ( 1/x4 - 1/x5 - p1.p1 );
        endif;
        #call ACCU(BN3 2)

        if( count(intbn3,1));
        id  p1.p6 = 1/2 * (-1/x3 + 1/x6 + p1.p1 );
        id  p2.p3 = 1/2 * ( 1/x5 - 1/x3 - p2.p2 );
        endif;
        #call ACCU(BN3 3)

        if( count(intbn3,1));
        id  p2.p4 = 1/2 * (-1/x6 + 1/x4 + p2.p2 );
        id  p2.p5 = 1/2 * (-1/x3 + 1/x5 + p2.p2 );
        endif;
        #call ACCU(BN3 4)

        if( count(intbn3,1));
        id  p2.p6 = 1/2 * ( 1/x4 - 1/x6 - p2.p2 );
        id  p3.p4 = 1/2 * ( 1/x5 + 1/x6 - p2.p2 - p1.p1 - 2*M^2);
        endif;
        #call ACCU(BN3 5)

        if( count(intbn3,1));
        id  p3.p5 = 1/2 * ( 1/x3 + 1/x5 - p2.p2 - 2*M^2);
        id  p3.p6 = 1/2 * ( 1/x3 + 1/x6 - p1.p1 - 2*M^2);
        id  p4.p5 = 1/2 * ( 1/x4 + 1/x5 - p1.p1 - 2*M^2);
        endif;
        #call ACCU(BN3 6)

        if( count(intbn3,1));
        id  p4.p6 = 1/2 * ( 1/x4 + 1/x6 - p2.p2 - 2*M^2);
        id  p5.p6 = 1/2 * ( 1/x3 + 1/x4 - p2.p2 - p1.p1 - 2*M^2);
        endif;
        #call ACCU(BN3 7)

*
* Warning! 
*
        if( count(intbn3,1));
        id  1/x3 = M^2 + p3.p3;
        id  1/x4 = M^2 + p4.p4;
        id  1/x5 = M^2 + p5.p5;

        id  p6.p6 = 1/x6 - M^2;
        endif;
        #call ACCU(BN3 5)

        #message do recursion

* no massive denominator:

        if( count(intbn3,1));
        if ( (count(x6,1)<=0) ) discard;
        endif;

        #call symBN3(p1,p2,p3,p4,p5,x6)

        if( count(intbn3,1));
        if ( (count(p5.p5,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        endif;
        .sort

* use the procedure triangle from MINCER to reduce the lines
* 2, 4 and 5.

* in order to avoid 1/fac_(-x)-terms:
        if( count(intbn3,1));
        if (count(p3.p3,1)>=0) multiply replace_(p1,p2,p2,p1,p4,p3,p3,p4);
        endif;
***if (count(p1.p1,1)>=0) multiply replace_(p1,p2,p2,p1,p4,p3,p3,p4);
        .sort

        if( count(intbn3,1));
        if (match(1/p2.p2/p4.p4/p5.p5)>0);
*         #call triown(p5,p3,p1,p2,p4)
        #call ntriangle(p5,p3,p1,p2,p4)
        endif;




        if ( (count(p5.p5,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        endif;
        .sort

*         Print+s;
*         .end

        if( count(intbn3,1));
        if ( ( match(1/p1.p1/p2.p2/p5.p5*x6) > 0 ) && (count(p4.p4,1)>=0) );
        id p4=p1+p5;
        if ( match(1/p1.p1/p2.p2/p5.p5*x6) <= 0 ) discard;
        multiply replace_(p1,p5,p5,p2,p2,p1);
        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6) > 0 ) && (count(p5.p5,1)>=0) );
        id p5=p4-p1;
        if ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6) <= 0 ) discard;
        id p1=-p1;
        id p4=-p4;
        multiply replace_(p3,p2,p2,p3);
        elseif ( ( match(1/p3.p3/p4.p4/p5.p5*x6) > 0 ) && (count(p2.p2,1)>=0) );
        id p2=p5-p3;
        if ( match(1/p3.p3/p4.p4/p5.p5*x6) <= 0 ) discard;
        multiply replace_(p3,p5,p1,p3,p4,p2,p5,p1);
        else;

* If more massless lines than in the three cases above are >= zero,
* a massless tadpole appears.

***  discard;
        multiply 1/(1-1);
        endif;

        endif;
        .sort

        #call Conv2exact

*         Print+s;
*         .end
        #message - done
#endprocedure



************************************************************

#procedure symBM (p1,p2,p3,x4,x5,x6)

        if( count(intbm,1));
        if (count(`x4',1) > count(`x6',1)) 
                    multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');
if (count(`x5',1) > count(`x6',1)) 
                    multiply replace_(`x5',`x6',`x6',`x5',`p1',`p2',`p2',`p1');
if (count(`x4',1) > count(`x5',1)) 
        multiply replace_(`x4',`x5',`x5',`x4',`p1',`p3',`p3',`p1');
endif;        
#endprocedure




#procedure symBMnom (p1,p2,p3,p4,p5,p6,x4,x5,x6)

* sort: n6>=n5>=n4
if( count(intbm,1));        
if ( (count(`x4',1) > count(`x6',1)) && (count(`x4',1) > count(`x5',1)) );
  id `p1'=-`p1';
  multiply replace_(`x4',`x6',`x6',`x4',
                    `p4',`p6',`p6',`p4',
                    `p2',`p3',`p3',`p2');
endif;
if ( (count(`x5',1) > count(`x6',1)) && (count(`x5',1) > count(`x4',1)) );
  id `p3'=-`p3';
  multiply replace_(`x5',`x6',`x6',`x5',
                    `p5',`p6',`p6',`p5',
                    `p1',`p2',`p2',`p1');
endif;

if (count(`x4',1) > count(`x5',1));
  id `p4'=-`p4';
  id `p5'=-`p5';
  id `p6'=-`p6';
  multiply replace_(`x5',`x4',`x4',`x5',
                    `p5',`p4',`p4',`p5',
                    `p3',`p1',`p1',`p3');
endif;
endif;
#endprocedure


#procedure nomBM
*
* Decomposition of the numerator for type BM
*

************************************************************

#call symBMnom(p1,p2,p3,p4,p5,p6,x4,x5,x6)
.sort

************************************************************

if( count(intbm,1));        
id p4.p4 = 1/x4 - M^2;
id p5.p5 = 1/x5 - M^2;
id p6.p6 = 1/x6 - M^2;

if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
  discard;
if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
endif;

#call ACCU(pi.pi)

************************************************************

if( count(intbm,1));        
id  p1.p2 = 1/2 * (-p3.p3 + p1.p1 + p2.p2);
id  p1.p3 = 1/2 * ( p2.p2 - p1.p1 - p3.p3);
id  p1.p4 = 1/2 * (-1/x6  + 1/x4  + p1.p1);
id  p1.p5 = 1/2 * ( p3.p3 + 1/x4  - p2.p2 - 1/x6);

if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
  discard;
if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
endif;

#call ACCU(p1)

************************************************************

if( count(intbm,1));        
id  p1.p6 = 1/2 * ( 1/x4  - p1.p1 - 1/x6);
id  p2.p3 = 1/2 * (-p1.p1 + p2.p2 + p3.p3);
id  p2.p4 = 1/2 * (-1/x5  + p2.p2 + 1/x4);
id  p2.p5 = 1/2 * ( 1/x4  - p2.p2 - 1/x5);

if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
  discard;
if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
endif;

#call ACCU(p1 p2)

************************************************************

if( count(intbm,1));        
id  p2.p6 = 1/2 * ( p3.p3 + 1/x4  - p1.p1 - 1/x5);
id  p3.p4 = 1/2 * ( p2.p2 + 1/x6  - p1.p1 - 1/x5);
id  p3.p5 = 1/2 * ( 1/x6  - p3.p3 - 1/x5);
id  p3.p6 = 1/2 * (-1/x5  + p3.p3 + 1/x6);

if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
  discard;
if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
endif;

#call ACCU(p2 p3)

************************************************************

if( count(intbm,1));        
id  p4.p5 = 1/2 * (-p2.p2 + 1/x4  + 1/x5 - 2*M^2);
id  p4.p6 = 1/2 * (-p1.p1 + 1/x4  + 1/x6 - 2*M^2);
id  p5.p6 = 1/2 * (-p3.p3 + 1/x5  + 1/x6 - 2*M^2);

if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
  discard;
if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
endif;

#call ACCU(p4 p5 p6)

************************************************************

#call symBM(p1,p2,p3,x4,x5,x6)

***#message numerator decomposition done (BM)
* #include matad.info # time
* #include matad.info # print

************************************************************
        
#endprocedure        






#procedure redBMn456 (p1,p2,p3,x4,x5,x6)

* Warning! In this procedure x3 is used for historical reasons!
*          Therefore not with `x3'

if( count(intbm,1));        
repeat;

* sort: n6>=n4,n5

  if((count(`x4',1) > 0)  &&  (count(x3,1) = 0)  &&
     (count(`x5',1) > 0)  &&  (count(`x6',1) > 0) );

    if (count(`x4',1) > count(`x6',1)) 
                    multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');
    if (count(`x5',1) > count(`x6',1)) 
                    multiply replace_(`x5',`x6',`x6',`x5',`p1',`p2',`p2',`p1');

* Now do the reduction via eq. (M7) resp. (4)

    if (count(`x6',1)>1);

      id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3?
        * `x4'^n4? * `x5'^n5? * `x6'^n6?
      =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3
        * `x4'^n4  * `x5'^n5  * `x6'^(n6-1) * (-1)
        *1/(n6-1)/M^2/2
        *(
          nom(4-2*(n6-1)-n1-n3,-2)
         +n1/`p1'.`p1'*(1/`x4'-1/`x6')
         +n3/`p3'.`p3'*(1/`x5'-1/`x6')
         )
        ;
    endif;

  endif;

endrepeat;
endif;

#endprocedure

#procedure redBMn123 (p1,p2,p3,x4,x5,x6)

* Warning! In this procedure x3 is used for historical reasons!
*          Therefore not with `x3'

* at this stage: n4=n5=n6=1,0
*                n1,n2,n3 <>=0

if( count(intbm,1));                
  repeat;

    if (count(`p3'.`p3',1) > count(`p1'.`p1',1)) 
                    multiply replace_(`x4',`x5',`x5',`x4',`p1',`p3',`p3',`p1');
    if (count(`p3'.`p3',1) > count(`p2'.`p2',1)) 
                    multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');


    if ( (count(`p3'.`p3',1) < 0) 
      && (count(`x4',1) = 1)       &&  (count(x3,1) = 0)  
      && (count(`x5',1) = 1)       &&  (count(`x6',1) = 1) );  

       id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
        * `x4'^n4? * `x5'^n5? * `x6'^n6?
      =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
        * `x4'^n4  * `x5'^n5  * `x6'^n6
        *1/2*deno(n1+n2-3,2)*deno(3-2*n3,-2)*deno(n1+n2+n3-3,2)
        *(-1)
        *(
        -(`x5'*nom(-3 + n1 + n2,2)**2/`x4') - `x6'*nom(-3 + n1 + n2,2)**2/`x4' 
+(-1)*
        n3*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x5'*`p3'.`p3') 
-(-1)*
        n3*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x6'*`p3'.`p3') 
-(-1)*
        n3*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x5'*`p3'.`p3') 
+(-1)*
        n3*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x6'*`p3'.`p3') 
+(-1)*
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(3 - n1 - n3,-2)/
         (2*M**2*`x5'*`p2'.`p2') 
+(-1)*
        n2*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(3 - n1 - n3,-2)/
         (2*M**2*`x4'*`p2'.`p2') 
+(-1)*
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(3 - n2 - n3,-2)/
         (2*M**2*`x6'*`p1'.`p1') 
+(-1)*
        n1*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(3 - n2 - n3,-2)/
         (2*M**2*`x4'*`p1'.`p1') 
+
        `x4'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/`x6' 
+
        `x5'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/`x6' 
+(-1)*
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/
         (2*M**2*`x4'*`p2'.`p2') 
+(-1)*
        n2*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/
         (2*M**2*`x5'*`p2'.`p2') 
+
        `x4'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/`x5' 
+
        `x6'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/`x5' 
+(-1)*
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/
         (2*M**2*`x4'*`p1'.`p1') 
+(-1)*
        n1*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/
         (2*M**2*`x6'*`p1'.`p1') 
+(-1)*
          `p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)*
          nom(-3 + n2 + n3,2)/M**2 
+(-1)*
        nom(-3 + n1 + n2,2)*nom(-6 + 9*n2 - n1*n2 - 2*n2**2 + n3 + n1*n3 -
            2*n2*n3,4 - 4*n2)/(2*M**2*`x5') 
+(-1)*
        nom(-3 + n1 + n2,2)*nom(-6 + 9*n1 - 2*n1**2 - n1*n2 + n3 - 2*n1*n3 +
            n2*n3,4 - 4*n1)/(2*M**2*`x6') 
+(-1)*
        nom(-3 + n1 + n2,2)*nom(12 - 9*n1 + 2*n1**2 - 9*n2 + 2*n1*n2 +
            2*n2**2 - 2*n3 + n1*n3 + n2*n3,-8 + 4*n1 + 4*n2)/(2*M**2*`x4')
      );
    endif;

    id nom(x?,0)=x;
    if ( (count(`x4',1)==0) && (count(x3,1)==0) && 
       ((count(`p2'.`p2',1)>=0) || (count(`p1'.`p1',1)>=0)) ) discard;

endrepeat;
endif;
  .sort

if( count(intbm,1));        
  repeat;

    if (count(`p3'.`p3',1) < count(`p1'.`p1',1)) 
                    multiply replace_(`x4',`x5',`x5',`x4',`p1',`p3',`p3',`p1');
    if (count(`p3'.`p3',1) < count(`p2'.`p2',1)) 
                    multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');

    if ( (count(`p3'.`p3',1) > 0)
      && (count(`x4',1) = 1)       &&  (count(x3,1) = 0)  
      && (count(`x5',1) = 1)       &&  (count(`x6',1) = 1) );  

       id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
        * `x4'^n4? * `x5'^n5? * `x6'^n6?
      =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
        * `x4'^n4  * `x5'^n5  * `x6'^n6
        *deno(n1+n2-3,2)*deno(n1+n3-2,2)*deno(n2+n3-2,2)
        *(-1)
        *(
       (1 + n3)*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*`x5'*`p3'.`p3'**2) 
-
        (1 + n3)*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*`x6'*`p3'.`p3'**2) 
-
        (1 + n3)*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*`x5'*`p3'.`p3'**2) 
+
        (1 + n3)*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*`x6'*`p3'.`p3'**2) 
+
        M**2*`x5'*nom(-3 + n1 + n2,2)**2/(`x4'*`p3'.`p3') 
+
        M**2*`x6'*nom(-3 + n1 + n2,2)**2/(`x4'*`p3'.`p3') 
+
        n2*nom(-3 + n1 + n2,2)*nom(2 - n1 - n3,-2)/(2*`x4'*`p2'.`p2') 
+
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(2 - n1 - n3,-2)/
         (2*`x5'*`p2'.`p2'*`p3'.`p3') 
+
        n1*nom(-3 + n1 + n2,2)*nom(2 - n2 - n3,-2)/(2*`x4'*`p1'.`p1') 
+
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(2 - n2 - n3,-2)/
         (2*`x6'*`p1'.`p1'*`p3'.`p3') 
+
        n2*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(2*`x5'*`p2'.`p2') 
-
        M**2*`x4'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(`x6'*`p3'.`p3') 
-
        M**2*`x5'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(`x6'*`p3'.`p3') 
+
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/
         (2*`x4'*`p2'.`p2'*`p3'.`p3') 
+
        n1*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(2*`x6'*`p1'.`p1') 
-
        M**2*`x4'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(`x5'*`p3'.`p3') 
-
        M**2*`x6'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(`x5'*`p3'.`p3') 
+
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/
         (2*`x4'*`p1'.`p1'*`p3'.`p3') 
-
        2*M**2*nom(-3 + n1 + n2,2)*nom(1 - 2*n3,-2)*nom(-2 + n1 + n2 + n3,2)/
         `p3'.`p3' 
+       nom(-3 + n1 + n2,2)*
          nom(-5 + n1 + 7*n2 - n1*n2 - 2*n2**2 + n3 + n1*n3 - 2*n2*n3,
           4 - 4*n2)/(2*`x5'*`p3'.`p3') 
+
        nom(-3 + n1 + n2,2)*nom(-5 + 7*n1 - 2*n1**2 + n2 - n1*n2 + n3 -
            2*n1*n3 + n2*n3,4 - 4*n1)/(2*`x6'*`p3'.`p3') 
+
        nom(-3 + n1 + n2,2)*nom(10 - 8*n1 + 2*n1**2 - 8*n2 + 2*n1*n2 +
            2*n2**2 - 2*n3 + n1*n3 + n2*n3,-8 + 4*n1 + 4*n2)/
         (2*`x4'*`p3'.`p3')
       );

    endif;

    id nom(x?,0)=x;
    if ( (count(`x4',1)==0) && (count(x3,1)==0) && 
       ((count(`p2'.`p2',1)>=0) || (count(`p1'.`p1',1)>=0)) ) discard;

  endrepeat;
  endif;
  
#endprocedure

#procedure redBMn25 (p1,p2,p3,x4,x5,x6)

if( count(intbm,1));                
repeat;
  if ( (count(`p2'.`p2',1) < 0) && (count(`p1'.`p1',1) == 0)
    &&  (count(`x4',1) >= 1)  
    && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  

     id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
      * `x4'^n4? * `x5'^n5? * `x6'^n6?
    =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
      * `x4'^n4  * `x5'^n5  * `x6'^n6
      *deno(4-2*n2-n4,-2)
      *(
        n4*`x4'*( `p2'.`p2' - 1/`x5' )
       )
      ;
  endif;

endrepeat;
endif;

#endprocedure

#procedure redBMn3 (p1,p2,p3,x4,x5,x6)

* Warning! In this procedure x3 is used for historical reasons!
*          Therefore not with `x3'

* at this stage: n4=n5=n6=1,0
*                n1,n2,n3 <>=0

if( count(intbm,1));                
repeat;

*if (count(`p3'.`p3',1) > count(`p1'.`p1',1)) 
*                 multiply replace_(`x4',`x5',`x5',`x4',`p1',`p3',`p3',`p1');
*if (count(`p3'.`p3',1) > count(`p2'.`p2',1)) 
*                 multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');

  if ( (count(`p3'.`p3',1) < 0)
     && (count(`p1'.`p1',1) != 0) &&  (count(`p2'.`p2',1) != 0)
     && (count(`x4',1) = 1)       &&  (count(x3,1) = 0)  
     && (count(`x5',1) = 1)       &&  (count(`x6',1) = 1) );  

     id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
      * `x4'^n4? * `x5'^n5? * `x6'^n6?
     =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
      * `x4'^n4  * `x5'^n5  * `x6'^n6
      *1/2*deno(n1+n2-3,2)*deno(3-2*n3,-2)*deno(n1+n2+n3-3,2)
      *(-1)
      *(
        -(`x5'*nom(-3 + n1 + n2,2)**2/`x4') - `x6'*nom(-3 + n1 + n2,2)**2/`x4' 
+(-1)*
        n3*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x5'*`p3'.`p3') 
-(-1)*
        n3*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x6'*`p3'.`p3') 
-(-1)*
        n3*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x5'*`p3'.`p3') 
+(-1)*
        n3*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*M**2*`x6'*`p3'.`p3') 
+(-1)*
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(3 - n1 - n3,-2)/
         (2*M**2*`x5'*`p2'.`p2') 
+(-1)*
        n2*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(3 - n1 - n3,-2)/
         (2*M**2*`x4'*`p2'.`p2') 
+(-1)*
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(3 - n2 - n3,-2)/
         (2*M**2*`x6'*`p1'.`p1') 
+(-1)*
        n1*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(3 - n2 - n3,-2)/
         (2*M**2*`x4'*`p1'.`p1') 
+
        `x4'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/`x6' 
+
        `x5'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/`x6' 
+(-1)*
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/
         (2*M**2*`x4'*`p2'.`p2') 
+(-1)*
        n2*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)/
         (2*M**2*`x5'*`p2'.`p2') 
+
        `x4'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/`x5' 
+
        `x6'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/`x5' 
+(-1)*
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/
         (2*M**2*`x4'*`p1'.`p1') 
+(-1)*
        n1*`p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n2 + n3,2)/
         (2*M**2*`x6'*`p1'.`p1') 
+(-1)*
          `p3'.`p3'*nom(-3 + n1 + n2,2)*nom(-3 + n1 + n3,2)*
          nom(-3 + n2 + n3,2)/M**2 
+(-1)*
        nom(-3 + n1 + n2,2)*nom(-6 + 9*n2 - n1*n2 - 2*n2**2 + n3 + n1*n3 -
            2*n2*n3,4 - 4*n2)/(2*M**2*`x5') 
+(-1)*
        nom(-3 + n1 + n2,2)*nom(-6 + 9*n1 - 2*n1**2 - n1*n2 + n3 - 2*n1*n3 +
            n2*n3,4 - 4*n1)/(2*M**2*`x6') 
+(-1)*
        nom(-3 + n1 + n2,2)*nom(12 - 9*n1 + 2*n1**2 - 9*n2 + 2*n1*n2 +
            2*n2**2 - 2*n3 + n1*n3 + n2*n3,-8 + 4*n1 + 4*n2)/(2*M**2*`x4')
      );
  endif;

  id nom(x?,0)=x;
  if ( (count(`x4',1)==0) && (count(x3,1)==0) && 
     ((count(`p2'.`p2',1)>=0) || (count(`p1'.`p1',1)>=0)) ) discard;

endrepeat;
endif;

#endprocedure

#procedure redBMn35 (p1,p2,p3,x4,x5,x6)

if( count(intbm,1));                
repeat;
  if ( (count(`p3'.`p3',1) < 0) && (count(`p1'.`p1',1) == 0)
    &&  (count(`x4',1) >= 1)  
    && (count(`x5',1) >= 1)  &&  (count(`x6',1) >= 1) );  

     id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
      * `x4'^n4? * `x5'^n5? * `x6'^n6?
     =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
      * `x4'^n4  * `x5'^n5  * `x6'^n6
      *deno(4-2*n3-n6,-2)
      *(
        n6*`x6'*( `p3'.`p3' - 1/`x5' )
       )
      ;
  endif;

endrepeat;
endif;

#endprocedure

#procedure redBMn3p (p1,p2,p3,x4,x5,x6)

if( count(intbm,1));                
repeat;

  if (count(`p3'.`p3',1) < count(`p1'.`p1',1)) 
                    multiply replace_(`x4',`x5',`x5',`x4',`p1',`p3',`p3',`p1');
  if (count(`p3'.`p3',1) < count(`p2'.`p2',1)) 
                    multiply replace_(`x4',`x6',`x6',`x4',`p2',`p3',`p3',`p2');

  if ( (count(`p3'.`p3',1) > 0)
    && (count(`x4',1) = 1)       &&  (count(x3,1) = 0)  
    && (count(`x5',1) = 1)       &&  (count(`x6',1) = 1) );  

     id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3? 
      * `x4'^n4? * `x5'^n5? * `x6'^n6?
     =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3 
      * `x4'^n4  * `x5'^n5  * `x6'^n6
      *deno(n1+n2-3,2)*deno(n1+n3-2,2)*deno(n2+n3-2,2)
      *(-1)
      *(
       (1 + n3)*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*`x5'*`p3'.`p3'**2) 
-
        (1 + n3)*`p1'.`p1'*nom(-3 + n1 + n2,2)**2/(2*`x6'*`p3'.`p3'**2) 
-
        (1 + n3)*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*`x5'*`p3'.`p3'**2) 
+
        (1 + n3)*`p2'.`p2'*nom(-3 + n1 + n2,2)**2/(2*`x6'*`p3'.`p3'**2) 
+
        M**2*`x5'*nom(-3 + n1 + n2,2)**2/(`x4'*`p3'.`p3') 
+
        M**2*`x6'*nom(-3 + n1 + n2,2)**2/(`x4'*`p3'.`p3') 
+
        n2*nom(-3 + n1 + n2,2)*nom(2 - n1 - n3,-2)/(2*`x4'*`p2'.`p2') 
+
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(2 - n1 - n3,-2)/
         (2*`x5'*`p2'.`p2'*`p3'.`p3') 
+
        n1*nom(-3 + n1 + n2,2)*nom(2 - n2 - n3,-2)/(2*`x4'*`p1'.`p1') 
+
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(2 - n2 - n3,-2)/
         (2*`x6'*`p1'.`p1'*`p3'.`p3') 
+
        n2*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(2*`x5'*`p2'.`p2') 
-
        M**2*`x4'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(`x6'*`p3'.`p3') 
-
        M**2*`x5'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/(`x6'*`p3'.`p3') 
+
        n2*`p1'.`p1'*nom(-3 + n1 + n2,2)*nom(-2 + n1 + n3,2)/
         (2*`x4'*`p2'.`p2'*`p3'.`p3') 
+
        n1*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(2*`x6'*`p1'.`p1') 
-
        M**2*`x4'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(`x5'*`p3'.`p3') 
-
        M**2*`x6'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/(`x5'*`p3'.`p3') 
+
        n1*`p2'.`p2'*nom(-3 + n1 + n2,2)*nom(-2 + n2 + n3,2)/
         (2*`x4'*`p1'.`p1'*`p3'.`p3') 
-
        2*M**2*nom(-3 + n1 + n2,2)*nom(1 - 2*n3,-2)*nom(-2 + n1 + n2 + n3,2)/
         `p3'.`p3' 
+       nom(-3 + n1 + n2,2)*
          nom(-5 + n1 + 7*n2 - n1*n2 - 2*n2**2 + n3 + n1*n3 - 2*n2*n3,
           4 - 4*n2)/(2*`x5'*`p3'.`p3') 
+
        nom(-3 + n1 + n2,2)*nom(-5 + 7*n1 - 2*n1**2 + n2 - n1*n2 + n3 -
            2*n1*n3 + n2*n3,4 - 4*n1)/(2*`x6'*`p3'.`p3') 
+
        nom(-3 + n1 + n2,2)*nom(10 - 8*n1 + 2*n1**2 - 8*n2 + 2*n1*n2 +
            2*n2**2 - 2*n3 + n1*n3 + n2*n3,-8 + 4*n1 + 4*n2)/
         (2*`x4'*`p3'.`p3')
       );

  endif;

  id nom(x?,0)=x;
  if ( (count(`x4',1)==0) && (count(x3,1)==0) && 
     ((count(`p2'.`p2',1)>=0) || (count(`p1'.`p1',1)>=0)) ) discard;

endrepeat;
endif;

#endprocedure

#procedure redBMn6 (p1,p2,p3,x4,x5,x6)

#do i = 1,1

* Now do the reduction via eq. (M7) resp. (4)

        if( count(intbm,1));        
  if ( (count(`x4',1)>=1) && (count(`x5',1)>=1) && (count(`x6',1)>1) );

    id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3?
      * `x4'^n4? * `x5'^n5? * `x6'^n6?
    =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3
      * `x4'^n4  * `x5'^n5  * `x6'^(n6-1) * (-1)
      *1/(n6-1)/M^2/2
      *(
        nom(4-2*(n6-1)-n1-n3, -2)
       +n1/`p1'.`p1'*(1/`x4'-1/`x6')
       +n3/`p3'.`p3'*(1/`x5'-1/`x6')
       )
      ;

    redefine i "0";

        endif;
endif;        
  .sort

#enddo

#endprocedure

#procedure redBMn6exp (p1,p2,p3,x4,x5,x6)

#do i = 1,1

* Now do the reduction via eq. (M7) resp. (4)

if( count(intbm,1));                
  if ( (count(`x4',1)>=1) && (count(`x5',1)>=1) && (count(`x6',1)>1) );

    id 1/`p1'.`p1'^n1? *1/`p2'.`p2'^n2? *1/`p3'.`p3'^n3?
      * `x4'^n4? * `x5'^n5? * `x6'^n6?
    =  1/`p1'.`p1'^n1 *1/`p2'.`p2'^n2 *1/`p3'.`p3'^n3
      * `x4'^n4  * `x5'^n5  * `x6'^(n6-1) * (-1)
      *1/(n6-1)/M^2/2
      *(
        num(4-2*(n6-1)-n1-n3 -2*ep)
       +n1/`p1'.`p1'*(1/`x4'-1/`x6')
       +n3/`p3'.`p3'*(1/`x5'-1/`x6')
       )
      ;
    redefine i "0";

  endif;
        endif;
        
*   repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
  #call ACCU{BMn6}

#enddo

#endprocedure


#procedure symmetryBM
*
* symmetryBM
*
if( count(intbm,1));                        
if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
  discard;
if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
  discard;
if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
endif;        
.sort

*
* sort n6>=n5>=n4
*
if( count(intbm,1));                        
if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
                             multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
                             multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) ) multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
endif;        
.sort

if( count(intbm,1));                        
        if ( (count(x4,1)<=0) || (count(x5,1)<=0) || (count(x6,1)<=0) );
endif;        
*
* sort n6>n5>n4
*
if( count(intbm,1));                        
  if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
                             multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
  if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
                             multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
  if ( count(x5,1) < count(x4,1) )
                             multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;
endif;        
.sort
        
#endprocedure        


#procedure reduceBMBM
*
* reduceBMBM
*

* This program reduces the integrals of type BM to easy integrals of
* the type M1 and T.

************************************************************


************************************************************

* massive lines

        #call symmetryBM
* #include expandnomdeno

* reduce n6 (expand in ep)
        #call redBMn6exp(p1,p2,p3,x4,x5,x6)
        .sort

        #call symmetryBM
* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

* now: n4=1 (or one of the massive lines = 0)

* #include expandnomdeno
        #call redBMn6exp(p1,p2,p3,x4,x5,x6)
        .sort

        #call symmetryBM
* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

* now: n4=n5=1

* #include expandnomdeno

        #call redBMn6exp(p1,p2,p3,x4,x5,x6)
        .sort

        #call symmetryBM
* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

************************************************************

* massless lines

* first treat all (!) massless indices which are <0

        #call redBMn3p(p1,p2,p3,x4,x5,x6)
        .sort

        #call symmetryBM
* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

        if( count(intbm,1));                        
        if ( (count(x4,1)>0) && (count(x5,1)>0) && (count(x6,1)>0) );
        
* sort: |n3|<=|n1|,|n2| ?
        
        
        if (count(p3.p3,1) < count(p1.p1,1)) 
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        if (count(p3.p3,1) < count(p2.p2,1)) 
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        endif;
        endif;
* the following procedure is only used, if all 3 indices are <0
        
        #call redBMn3(p1,p2,p3,x4,x5,x6)
        .sort

        #call symmetryBM

        if( count(intbm,1));                                
        if ( (count(x4,1)>0) && (count(x5,1)>0) && (count(x6,1)>0) );
        if (count(p3.p3,1)==0) multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;
        endif;        
        .sort

* now one massless index is =0; here n1=0

* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

* Now treat the two other massless lines with simpler
* rec. relations.

        #call redBMn35(p1,p2,p3,x4,x5,x6)
        .sort

* Don't use the file 'symmetryBM' here because the next procedure 
* treats n2>0.

        if( count(intbm,1));                                
        if ( (count(x4,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p2.p2,-1)<=0) ) )
        discard;
        if ( (count(x5,1)<=0) && ( (count(p2.p2,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x6,1)<=0) && ( (count(p1.p1,-1)<=0) || (count(p3.p3,-1)<=0) ) )
        discard;
        if ( (count(x4,1)<=0) && (count(x5,1)<=0) ) discard;
        if ( (count(x4,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        endif;        
        .sort

        if( count(intbm,1));                                
        if ( (count(x4,1)<=0) || (count(x5,1)<=0) || (count(x6,1)<=0) );
*
* sort n6>n5>n4
*
        if ( (count(x4,1) > count(x6,1)) && (count(x4,1) >= count(x5,1)) )
        multiply replace_(x4,x6,x6,x4,p2,p3,p3,p2);
        if ( (count(x5,1) > count(x6,1)) && (count(x5,1) >= count(x4,1)) )
        multiply replace_(x5,x6,x6,x5,p2,p1,p1,p2);
        if ( count(x5,1) < count(x4,1) )
        multiply replace_(x4,x5,x5,x4,p1,p3,p3,p1);
        endif;
        endif;        
        .sort

* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

        #call redBMn25(p1,p2,p3,x4,x5,x6)
        .sort

        #call symmetryBM
* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

* if n1=n2=n3=0 the following rec. rel. is very simple

        #CALL redBMn456(p1,p2,p3,x4,x5,x6)

        #call symmetryBM
* #include redcutnomdeno
* #include expandnomdeno

* #include matad.info # time
* #include matad.info # print

#endprocedure        



#procedure topbm
*
* this is topbmbm
*

************************************************************

        #message this is topbmbm
        
        #message numerator
        
        #call nomBM
        
        #message do recursion
        
        #call reduceBMBM
        
        if( count(intbm,1));        
        id x4^n4?neg_=(p4.p4+M^2)^-n4;
        id p4 = p2+p5;
        id p5.p5 = 1/x5 - M^2;
        
        if ( (count(x4,1)=0) && (count(x3,1)=0) 
        && ((count(p2.p2,1)>=0) || (count(p1.p1,1)>=0)) ) discard;
* this is necessary, because it is possible, that n1 and n2 shrink to
* a point!
        if (count(x5,1)<=0) discard;
        if (count(x6,1)<=0) discard;
        endif;        
        .sort

        #call Conv2exact()        

Print+s;
.end

        .sort

        #message - done
        
#endprocedure        



#procedure symBM1 (p1,p2,p3,p4,x5,x6)
        if( count(intbm1,1));
        if (count(`x6',1) > count(`x5',1))
        multiply replace_(`p1',`p2',`p2',`p1',`x5',`x6',`x6',`x5');
        endif;
#endprocedure


#procedure redBM1n12 (p1,p2,p3,p4,x5,x6)
        
        #do i=1,1
                if( count(intbm1,1));
                if ( count(`p1'.`p1',1) > count(`p2'.`p2',1) ) 
                multiply replace_(`p1',`p2',`p2',`p1',`x5',`x6',`x6',`x5');
                endif;
                #call ACCU{BM1}
                
                if( count(intbm1,1));        
                if (match(1/`p1'.`p1'^2/`p2'.`p2'/`p4'.`p4'*`x5'*`x6')>0);
                
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4  * `x5'^n5  * `x6'^n6
                *(-1)/(n1-1)/M^2*(
                - `p1'.`p1'*n3/`p3'.`p3'*(1/`x6'-1/`x5')
                - (n1-1)                *(1/`x6'-`p4'.`p4')
                + `p1'.`p1'*(n-2*n6-n3-(n1-1))
                +2*n6*M^2*`x6'*`p1'.`p1'
                );
                
                redefine i "0";
                
                endif;
                
                if (match(1/`p1'.`p1'/`p2'.`p2'^2/`p4'.`p4'*`x5'*`x6')>0);
                
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4  * `x5'^n5  * `x6'^n6
                *(-1)/(n2-1)/M^2*(
                - `p2'.`p2'*n3/`p3'.`p3'*(1/`x5'-1/`x6')
                - (n2-1)                *(1/`x5'-`p4'.`p4')
                + `p2'.`p2'*(n-2*n5-n3-(n2-1))
                +2*n5*M^2*`x5'*`p2'.`p2'
                );
                
                redefine i "0";

                endif;

                id n = num(4-2*ep);

* topBM1        
                endif;
                #call ACCU{BM1}
                
        #enddo
        
* #include expandnomdeno
        
#endprocedure



#procedure redBM1n124 (p1,p2,p3,p4,x5,x6)

        #do i=1,1
                if( count(intbm1,1));
                if (match(1/`p1'.`p1'/`p2'.`p2'/`p4'.`p4'*`x5'*`x6')>0);
                
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4  * `x5'^n5  * `x6'^n6
                *(-1)*deno(2+n4-n1-n2-n3,-1)
                *(
                n5*`x5'*(`p4'.`p4'-`p2'.`p2')
                +n6*`x6'*(`p4'.`p4'-`p1'.`p1')
                )
                ;
                
                redefine i "0";
                
                endif;
                endif;
                #call ACCU{BM1}
                
        #enddo
        
* #include expandnomdeno

#endprocedure



#procedure redBM1n3 (p1,p2,p3,p4,x5,x6)

        if( count(intbm1,1));
        if ( (match(1/`p2'.`p2'/`p3'.`p3'/`p4'.`p4'*`x5'*`x6')>0) 
        && (count(`p1'.`p1',1)>=0) );
        if ( (count(`p3'.`p3',1) < -1 ) && (count(`p1'.`p1',1) > 0 ) );
        id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
        /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
        =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
        /`p3'.`p3'^n3  /`p4'.`p4'^n4 * `x5'^n5 * `x6'^n6
        *(-1)/(n3-1)
        *(
        `p3'.`p3'/`p1'.`p1'*nom(4-2*n2-(n3-1)-n5,-2)
        -(n3-1)*`p2'.`p2'/`p1'.`p1'
        +`p3'.`p3'/`p1'.`p1'*n5*`x5'*(`p4'.`p4'-`p2'.`p2'+M^2)
        )
        ;
        
        redefine i "0";
        
        endif;
        endif;
        
        if ( (match(1/`p1'.`p1'/`p3'.`p3'/`p4'.`p4'*`x5'*`x6')>0) 
        && (count(`p2'.`p2',1)>=0) );
        if ( (count(`p3'.`p3',1) < -1 ) && (count(`p2'.`p2',1) > 0 ) );
        id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
        /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
        =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
        /`p3'.`p3'^n3  /`p4'.`p4'^n4 * `x5'^n5 * `x6'^n6
        *(-1)/(n3-1)
        *(
        `p3'.`p3'/`p2'.`p2'*nom(4-2*n1-(n3-1)-n6,-2)
        -(n3-1)*`p1'.`p1'/`p2'.`p2'
        +`p3'.`p3'/`p2'.`p2'*n6*`x6'*(`p4'.`p4'-`p1'.`p1'+M^2)
        )
        ;
        endif;
        endif;
* topBM1        
        endif;        
#endprocedure

#procedure redBM1n35 (p1,p2,p3,p4,x5,x6)
        
        #do i=1,10
                
                if( count(intbm1,1));        
                if ( (match(1/`p2'.`p2'/`p3'.`p3'/`p4'.`p4'*`x5'*`x6')>0) 
                && (count(`p1'.`p1',1)>=0) );
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4 * `x5'^n5 * `x6'^n6
                *(-1)*deno(4-2*n3-n1-n6,-2)
                *(
                n1/`p1'.`p1'*(`p2'.`p2'-`p3'.`p3')
                +n6*`x6'*(1/`x5'-`p3'.`p3')
                )
                ;
                
                redefine i "0";
                
                endif;
                endif;
                #call ACCU{BM1}
                
        #enddo
        
* #include expandnomdeno
        
#endprocedure


#procedure redBM1n36 (p1,p2,p3,p4,x5,x6)
        
* with this procedure n3 or n6 is reduced to zero.
        
        #do i=1,1
                
                if( count(intbm1,1));        
                if ( (match(1/`p1'.`p1'/`p3'.`p3'/`p4'.`p4'*`x5'*`x6')>0) 
                && (count(p2.p2,1)>=0) );
                
                id 1/`p1'.`p1'^n1? /`p2'.`p2'^n2? 
                /`p3'.`p3'^n3? /`p4'.`p4'^n4? * `x5'^n5? * `x6'^n6?
                =  1/`p1'.`p1'^n1  /`p2'.`p2'^n2 
                /`p3'.`p3'^n3  /`p4'.`p4'^n4 * `x5'^n5 * `x6'^n6
                *(-1)*deno(4-2*n3-n2-n5,-2)
                *(
                n2/`p2'.`p2'*(`p1'.`p1'-`p3'.`p3')
                +n5*`x5'*(1/`x6'-`p3'.`p3')
                )
                ;
                
                redefine i "0";
                
                endif;
                endif;
                
                #call ACCU{BM1}
                
        #enddo
        
* #include expandnomdeno
        
#endprocedure



#procedure topbm1
*
* this is topbm1
*
        #-
        #message this is topbm1
        .sort
        
************************************************************
        
        
************************************************************
        
        #message numerator
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        ID  p6.p6 = 1/x6 - M^2;
        ID  p5.p5 = 1/x5 - M^2;
        endif;        
        #call ACCU(BM1 0)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        id p1=p2-p3;
        id p6=p3+p5;
        endif;        
        #call ACCU(BM1 1)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        id  p2.p3 = 1/2 * (-p1.p1 + p2.p2 + p3.p3);
        id  p2.p4 = 1/2 * (-1/x5  + p2.p2 + M^2+p4.p4);
        endif;        
        #call ACCU(BM1 2)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        id  p3.p4 = 1/2 * ( p2.p2 + 1/x6  - p1.p1 - 1/x5);
        endif;        
        #call ACCU(BM1 3)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm1,1));        
        id  p2.p5 = 1/2 * ( M^2+p4.p4  - p2.p2 - 1/x5);
        id  p3.p5 = 1/2 * ( 1/x6  - p3.p3 - 1/x5);
        endif;        
        #call ACCU(BM1 5)
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort
        
        
        #call ACCU(BM1 6)
        
        if( count(intbm1,1));
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;
        .sort
        if( count(intbm1,1));
        id  p4.p5 = 1/2 * (-p2.p2 + p4.p4 + 1/x5 - M^2);
        endif;        
        #call ACCU(BM1 7)

        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort

*
* Warning!
*
        if( count(intbm1,1));        
        ID  p6.p6 = 1/x6 - M^2;
        ID  p5.p5 = 1/x5 - M^2;
        endif;        
        #call ACCU(BM1 8)
        
        if( count(intbm1,1));       
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort

        #call ACCU(BM1)

        #message do recursion

        #call symBM1(p1,p2,p3,p4,x5,x6)
        #call redBM1n12(p1,p2,p3,p4,x5,x6)

        #call symBM1(p1,p2,p3,p4,x5,x6)
        #call redBM1n124(p1,p2,p3,p4,x5,x6)
        .sort

* from now on n1, n2 or n4 is <= zero. 
        if( count(intbm1,1));
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort

* #include expandnomdeno

        #do i=1,1
***  #message BM1-n3 'i'
                #call symBM1(p1,p2,p3,p4,x5,x6)
                #call redBM1n3(p1,p2,p3,p4,x5,x6)
                if( count(intbm1,1));        
                if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
                if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
                if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
                if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
                if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
                if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
                if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
                if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
                if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
                if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
                endif;        
                #call ACCU(BM1_n3)
        #enddo
        
        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        #call ACCU(BM1_n3)

        #call redBM1n35(p1,p2,p3,p4,x5,x6)

* #include expandnomdeno

        #call redBM1n36(p1,p2,p3,p4,x5,x6)

* #include expandnomdeno

        if( count(intbm1,1));        
        if ( (match(1/p1.p1/p2.p2*x5*x6)>0) && (count(p4.p4,1)>=0) );
        id p4=p2+p5;
        if (match(1/p1.p1/p2.p2*x5*x6)<=0) discard;
        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x5) > 0 )  && (count(x6,1)<=0) );
        id 1/x6 = M^2 + p6.p6;
        id p6=p3+p5;
        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6) > 0 )  && (count(x5,1)<=0) );
        id 1/x5 = M^2 + p5.p5;
        id p5=p2+p4;
***else;
***  multiply 1/(1-1);
        endif;
        endif;        
        .sort

* #include expandnomdeno

* treat negative powers of massive denominators:

        if( count(intbm1,1));        
        id x5^n5?neg_=(p4.p4-2*p4.p2+p2.p2+M^2)^-n5;
        id x6^n6?neg_=(p4.p4-2*p4.p1+p1.p1+M^2)^-n6;
        endif;        
        .sort

        if( count(intbm1,1));        
        if ( (count(x5,1)<=0) && (count(x6,1)<=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x5,1)<=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(x6,1)<=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        endif;        
        .sort

* now, we have reduced diagrams with the following lines != 0:

* 1,2,(3),5,6 : M1
* 2,4,5,6   : M4
* 2,3,4,6   : M2
* 1,4,5,6   : M4
* 1,3,4,5   : M2
* 1,2,3,4,5   : M2
* 1,2,3,4,6   : M2

* the notation has to be adjusted

        if( count(intbm1,1));        
        if ( (match(1/p2.p2/p4.p4*x5*x6) > 0) && (count(p1.p1,1)>=0) 
        && (count(p3.p3,1)>=0) );
        id p1=p4-p6;
        id p3=p6-p5;
        id p2=-p2;
        if ( match(1/p2.p2/p4.p4*x5*x6) <= 0 ) discard;
        multiply replace_(p6,p5,p5,p6,p2,p3,x6,x5,x5,x6);
        Multiply intm4/intbm1;

        elseif ( (match(1/p2.p2/p3.p3/p4.p4*x6) > 0)  && (count(p1.p1,1)>=0) 
        && (count(x5,1)<=0) );
        id p1=p4-p6;
        id p5=p6-p3;
        if ( match(1/p2.p2/p3.p3/p4.p4*x6) <= 0 ) discard;
        multiply replace_(p2,p1,p3,p2,p4,p5);
        Multiply intm2/intbm1;

        elseif ( (match(1/p1.p1/p4.p4*x5*x6) > 0)  && (count(p3.p3,1)>=0) 
        && (count(p2.p2,1)>=0) );
        id p2=p4-p5;
        id p3=p6-p5;
        id p6=-p6;
        id p4=-p4;
        if ( match(1/p1.p1/p4.p4*x5*x6) <= 0 ) discard;
        multiply replace_(p1,p3);
        Multiply intm4/intbm1;
        
        elseif ( (match(1/p1.p1/p3.p3/p4.p4*x5) > 0 ) && (count(p2.p2,1)>=0) 
        && (count(x6,1)<=0) );;
        id p2=p1+p3;
        id p6=p3+p5;
        id p3=-p3;
        if ( match(1/p1.p1/p3.p3/p4.p4*x5) <= 0 ) discard;
        multiply replace_(p5,p6,p3,p2,p4,p5,x5,x6);
        Multiply intm2/intbm1;

        elseif ( ( match(1/p1.p1/p2.p2*x5*x6) > 0 )  && (count(p4.p4,1)>=0) );

* This is type M1 and the notation is already o.k.
        Multiply intm1/intbm1;

        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x5) > 0 )  && (count(x6,1)<=0) );
        
* This is type M2

        id p1=-p1;
        id p2=-p2;
        id p3=-p3;
        id p4=-p4;
        id p5=-p5;
        multiply replace_(p4,p5,p1,p2,p3,p1,p2,p3,p5,p6,x5,x6);
        Multiply intm2/intbm1;

        elseif ( ( match(1/p1.p1/p2.p2/p3.p3/p4.p4*x6) > 0 )  && (count(x5,1)<=0) );

* This is type M2

        id p2=-p2;
        id p4=-p4;
        id p6=-p6;
        multiply replace_(p4,p5,p1,p3,p3,p1);
        Multiply intm2/intbm1;

        else;
        multiply 1/(1-1);
        endif;
* topBM1        
        endif;
        #call ACCU(BM1)

        #message - done
        
#endprocedure        


#procedure symBM2 (p1,p2,p3,p4,p5,x6)

        if( count(intbm2,1));        
        if (count(`p3'.`p3',1) >= count(`p1'.`p1',1))
        multiply replace_(`p1',`p3',`p3',`p1',`p4',`p5',`p5',`p4');
        endif;
#endprocedure


#procedure topbm2
*
* this is topbm2 
*
#-

#message this is topbm2

************************************************************


************************************************************

        #message numerator
        
        if( count(intbm2,1));        
        id  p6.p6 = 1/x6 - M^2;
        endif;        
        #call ACCU(BM2 0)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;
        .sort
        
        if( count(intbm2,1));                
        id  p1.p2 = 1/2 * (-p3.p3 + p1.p1 + p2.p2);
        id  p1.p3 = 1/2 * ( p2.p2 - p1.p1 - p3.p3);
        id  p1.p4 = 1/2 * (-1/x6  + 1/x4  + p1.p1);
        endif;        
        #call ACCU(BM2 1)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm2,1));                
        id  p1.p5 = 1/2 * ( p3.p3 + 1/x4  - p2.p2 - 1/x6);
        id  p1.p6 = 1/2 * ( 1/x4  - p1.p1 - 1/x6);
        id  p2.p3 = 1/2 * (-p1.p1 + p2.p2 + p3.p3);
        endif;        
        #call ACCU(BM2 2)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm2,1));                
        id  p2.p4 = 1/2 * (-1/x5  + p2.p2 + 1/x4);
        id  p2.p5 = 1/2 * ( 1/x4  - p2.p2 - 1/x5);
        id  p2.p6 = 1/2 * ( p3.p3 + 1/x4  - p1.p1 - 1/x5);
        endif;        
        #call ACCU(BM2 3)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort        
        
        if( count(intbm2,1));                
        id  p3.p4 = 1/2 * ( p2.p2 + 1/x6  - p1.p1 - 1/x5);
        id  p3.p5 = 1/2 * ( 1/x6  - p3.p3 - 1/x5);
        id  p3.p6 = 1/2 * (-1/x5  + p3.p3 + 1/x6);
        endif;        
        #call ACCU(BM2 4)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm2,1));                
        id  p4.p5 = 1/2 * (-p2.p2 + 1/x4  + 1/x5 - 2*M^2);
        id  p4.p6 = 1/2 * (-p1.p1 + 1/x4  + 1/x6 - 2*M^2);
        id  p5.p6 = 1/2 * (-p3.p3 + 1/x5  + 1/x6 - 2*M^2);
        endif;        
        #call ACCU(BM2 5)
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort

*
* Warning! 
*
        if( count(intbm2,1));                
        id  1/x4 = M^2 + p4.p4;
        id  1/x5 = M^2 + p5.p5;
        endif;        
        #call ACCU(BM2 6)
        
        if( count(intbm2,1));                
        id  p6.p6 = 1/x6 - M^2;
        endif;        
        #call ACCU(BM2 7)
***.SORT
        
        #message do recursion

* no massive denominator:
        
        if( count(intbm2,1));                
        if ( (count(x6,1)<=0) ) discard;
        endif;
        
        #call symBM2(p1,p2,p3,p4,p5,x6)
        
        if( count(intbm2,1));                
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
* Use the procedure triangle from MINCER to reduce the lines
* 1, 2 and 4.
        
* in order to avoid 1/fac_(-x)-terms:

*         Print+s;
*         .end        

        if( count(intbm2,1));                
        if (count(p5.p5,1)>=0) multiply replace_(p5,p4,p4,p5,p1,p3,p3,p1);
        
        if (match(1/p1.p1/p2.p2/p4.p4)>0);
***#call triangle(p2,p3,p5,p1,p4)
        #call triown(p2,p3,p5,p1,p4)
        endif;
        
        if ( (count(p5.p5,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p5.p5,1)>=0) && (count(p4.p4,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p1.p1,1)>=0) ) discard;
        if ( (count(p4.p4,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        if ( (count(p1.p1,1)>=0) && (count(p2.p2,1)>=0) ) discard;
        if ( (count(p2.p2,1)>=0) && (count(p3.p3,1)>=0) ) discard;
        endif;        
        .sort
        
        if( count(intbm2,1));        
        if ( ( match(1/p2.p2/p3.p3/p4.p4*x6) > 0 ) && (count(p1.p1,1)>=0) );
        id p1=p2-p3;
        if ( match(1/p2.p2/p3.p3/p4.p4*x6) <= 0 ) discard;
        multiply replace_(p5,p3,p3,p5,p4,p2,p2,p1);
        elseif ( ( match(1/p1.p1/p3.p3/p4.p4/p5.p5*x6) > 0 ) && (count(p2.p2,1)>=0) );
        id p2=p1+p3;
        id p5=-p5;
        if ( match(1/p1.p1/p3.p3/p4.p4/p5.p5*x6) <= 0 ) discard;
        multiply replace_(p4,p2,p5,p4);
        elseif ( ( match(1/p1.p1/p2.p2/p5.p5*x6) > 0 ) && (count(p4.p4,1)>=0) );
        id p4=p2+p5;
        if ( match(1/p1.p1/p2.p2/p5.p5*x6) <= 0 ) discard;
        else;
        
* If more massless lines than in the three cases above are >= zero,
* a massless tadpole appears.
        
***  discard;
        multiply 1/(1-1);
        endif;
        endif;        
        .sort
        
        #message - done
        
#endprocedure        



#procedure nomgm3(v1,v2,v3,x,y,z)

id  `v3' = -`v1'-`v2';

#call ACCU{nomgm3 1}

id  `v1'.`v1' = 1/`x' - M^2;

#call ACCU{nomgm3 2}

id  `v2'.`v2' = 1/`y' - M^2;

#call ACCU{nomgm3 3}

id  `v1'.`v2' = (1/(`z') - 1/`x' - 1/`y' +2*M^2)/2;

#call ACCU{nomgm3 4}

#endprocedure



#procedure topm1
*
* topm1
* (simple topology)
*

#message this is topm1
        
        if( count(intm1,1));                
*         multiply int1;
        id p1=-p1;
        endif;                
*         #call one00(p2,p1,p3,e2,e1,e3,int1)
        #call IntOne(p2,p1,p3,m1,MM0)
        .sort
        
#call ACCU(one)

* #call simpfin
* #call ACCU(simp)
* #include redcut

* #include expandnomdeno

        if( count(intMM0,1));
        id p6=-p6;
        endif;                
        #call nomgm3(p5,p6,p3,x5,x6,1/p3.p3)

        .sort

        Print+s;
        .sort: Before Tad2l;
        
*         #call two110(x5,s5m1,x6,s6m1,1/p3.p3,e3)
        #call TadpoleMM0(x5,x6,p3,MM0,0)
        .sort
        
        #call ACCU(gm3)
* #include expandnomdeno
        
#endprocedure        



#procedure topm2
*
* topm2
* (simple topology)
*
        
        #message this is topm2
        
        if( count(intm2,1));        
*         multiply int1;
        id p1=-p1;
        endif;        
*         #call one00(p2,p1,p3,e2,e1,e3,int1)
        #call IntOne(p2,p1,p3,m2,M00)        
        .sort

* #include expandnomdeno
        
        #call IntOne(p3,p5,p6,M00,M0)        
* multiply int1;
* #call one00(p3,p5,p6,e3,e5,e6,int1)
        .sort
        #call averts(p6,M0)
        .sort
*         #call one10(x6,s6m1,1/p6.p6,e6)
        #call TadpoleM0(x6,p6,M0,0)        
        .sort
        
* #include expandnomdeno
        
#endprocedure        




#procedure topm3
*
* topm3
* (simple topology)
*

        #message this is topm3

* multiply int1;
        if( count(intm3,1)) id p1=-p1;
*         #call one00(p2,p1,p6,e2,e1,e6,int1)
        #call IntOne(p2,p1,p6,m3,M00)                
.sort

* #include expandnomdeno

if( count(intm3,1)) id p4=-p4;
* multiply int1;
*         #call one00(p3,p4,p6,e3,e4,e6,int1)
        #call IntOne(p3,p4,p6,M00,M0)        
        .sort
#call averts(p6,M0)
.sort
*         #call one10(x6,s6m1,1/p6.p6,e6)
#call TadpoleM0(x6,p6,M0,0)        
.sort

* #include expandnomdeno
        
#endprocedure        



#procedure topm4
*
* topm4
* (simple topology)
*

#message this is topm4
* if( count(intm4,1)) multiply intM00*intM0/intm4;
*         #call one00(p3,p4,p6,e3,e4,e6,int1)
        #call IntOne(p3,p4,p6,m4,MxM)
.sort

* #include expandnomdeno

#call averts(p5,MxM)
.sort
*         #call one10(x5,s5m1,1/p5.p5,e5)
        #call TadpoleM0(x5,p5,MxM,M0)        
.sort
#call averts(p6,M0)
.sort
*         #call one10(x6,p6)
#call TadpoleM0(x6,p6,M0,0)        
.sort

* #include expandnomdeno
        
#endprocedure        


#procedure topm5
*
* this is topm5
*

        #message this is topm5

************************************************************

* treat 1-loop integral

        if( count(intm5,1)) multiply replace_(x1,s1m,x2,s2m,x3,s3m,x4,s4m,x5,s5m,x6,s6m,x7,s7m,x8,s8m);

        #call averts(p4,m5)
*         #call one10(s4m,s4m1,1/p4.p4,yy1)
        #call TadpoleM0(s4m,p4,m5,MMM)
        .sort

* #call simpfin
* #include redcut
* #include expandnomdeno

* treat numerator of 2-loop integral

        
        if( count(intMMM,1)) id p1.p1 = 1/s1m - M^2;
        .sort
        if( count(intMMM,1)) id p1 = p2+p5;
        .sort

        if( count(intMMM,1));        
        id p2.p2 = 1/s2m - M^2;
        id p5.p5 = 1/s5m - M^2;
        endif;        
        .sort

        if( count(intMMM,1)) id p2.p5 = 1/2 * (p1.p1-p2.p2-p5.p5);
        .sort

        if( count(intMMM,1));        
        id p1.p1 = 1/s1m - M^2;
        id p2.p2 = 1/s2m - M^2;
        id p5.p5 = 1/s5m - M^2;
        endif;        
        .sort

* recursion for 2-loop integral
        
* the following part is identical to topT1
        
        #do i=1,5
                
***#message 'i'
                if( count(intMMM,1));        
                if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s5m,1)!=0) ); 
                multiply replace_(test5,s1m,test6,s2m);
                if ( (count(s1m,1)>1) && (count(s2m,1)>=1) && (count(s5m,1)>=1) );
                id s1m^n1?*s2m^n2?*s5m^n5? = 
                s1m^n1*s2m^n2*s5m^n5 * (-1) * 1/3/(n1-1)/M^2 * (
                nom(4+3-3*n1,-2)/s1m
                +2*n2*s2m/s1m*(1/s5m-1/s1m)
                -(n1-1)*(1/s5m-1/s2m)
                );
                endif;
                if ( count(s1m,1) < count(s2m,1) ) multiply replace_(s1m,s2m,s2m,s1m);
                if ( count(s1m,1) < count(s5m,1) ) multiply replace_(s1m,s5m,s5m,s1m);
                if ( count(s2m,1) < count(s5m,1) ) multiply replace_(s2m,s5m,s5m,s2m);
                endif;
                endif;        
                .sort
                
* #include expandnomdeno
                
        #enddo

        if( count(intMMM,1));
        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s5m,1)!=0) ); 
        multiply replace_(test5,s1m,test6,s2m);
        repeat;
                if ( (count(s1m,1)>1) && (count(s2m,1)>=1) && (count(s5m,1)>=1) );
                id s1m^n1?*s2m^n2?*s5m^n5? = 
                s1m^n1*s2m^n2*s5m^n5 * (-1) * 1/3/(n1-1)/M^2 * (
                nom(4+3-3*n1,-2)/s1m
                +2*n2*s2m/s1m*(1/s5m-1/s1m)
                -(n1-1)*(1/s5m-1/s2m)
                );
                endif;
                if ( count(s1m,1) < count(s2m,1) ) multiply replace_(s1m,s2m,s2m,s1m);
                if ( count(s1m,1) < count(s5m,1) ) multiply replace_(s1m,s5m,s5m,s1m);
                if ( count(s2m,1) < count(s5m,1) ) multiply replace_(s2m,s5m,s5m,s2m);
        endrepeat;
        endif;
        endif;

        .sort

        #message perform integration

        if( count(intMMM,1));  
        if ( count(s1m,1) < count(s2m,1) ) multiply replace_(s1m,s2m,s2m,s1m);
        if ( count(s1m,1) < count(s5m,1) ) multiply replace_(s1m,s5m,s5m,s1m);
        if ( count(s2m,1) < count(s5m,1) ) multiply replace_(s2m,s5m,s5m,s2m);

        id 1/s5m = p5.p5 + M^2;
        endif;  
        .sort

        if( count(intMMM,1));  
        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s5m,1)==0) ); 
***  #call nomgm3T(p1,p2,p5,test5,test6,1/p5.p5)
        multiply replace_(test5,s1m,test6,s2m)*intMM0/intMMM;

        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s5m,1)==0) );
        multiply intM00/intm5;
        elseif ( (count(s1m,1)==1) && (count(s2m,1)==1) && (count(s5m,1)==1) ); 
***  id s1m*s2m*s5m = -3*deno(1,-2)*M^2*(1/2/ep^2 + 1/2/ep) +  M^2*g3fl1;
*** constant form Davydychev, Tausk, NPB397(93)123
        id s1m*s2m*s5m = M^2*(
        -21/2 - 3/(2*ep^2) - 9/(2*ep) + (27*S2)/2 - (3*z2)/2 
        + ep * T1ep + ep^2 * T1ep2
        );
        Multiply int0/intMMM;
        else;
        multiply 1/(1-1);
        endif;
        endif;
* MM0  
*   #call two110(s1m,s1m1,s2m,s2m1,1/p5.p5,e5)
        #call TadpoleMM0(s1m,s2m,p5,MM0,0)
* M00
*    multiply int1;
*   #call one00(p2,p5,p1,e2,e5,e1,int1)
        #call IntOne(p2,p5,p1,M00,M0)  
        #call averts(p1,M0)
*   #call one10(s1m,s1m1,1/p1.p1,e1)
        #call TadpoleM0(s1m,p1,M0,0)
        .sort

************************************************************
#endprocedure


#procedure topt1
#message this is topt1

        if( count(intt1,1));
        if( (count(x3,1)=0) && (count(x4,1)=1)
        && (count(x5,1)=1) && (count(x6,1)=1) )
        id x4*x5*x6 =  + (M^2*Gam(-1,1))^3;

*         Multiply int0/intt1;        
        endif;
     
*         #call GammaArgToOne   
* .sort
*         #call expansion(3)
*         b ep;
*         Print+s;
*         .end        

#endprocedure        




        
#procedure topn1
#message this is topn1

* Print+s;
* .end
************************************************************

* integrals of type N1
        
        if(count(intn1,1));        
        if( (count(x3,1)=1) && (count(x4,1)=1)
        && (count(x5,1)=1) && (count(x6,1)=1) )
        id x3*x4*x5*x6 =M^4*rat(agam/ep^3+bgam/ep^2+cgam/ep+dgam+egam*ep
        +fgam*ep^2+ggam*ep^3+hgam*ep^4,1);
        
*         Multiply int0/intn1;        
        endif;        
.sort 

* The following commands should not be used if 'EXACTEXP' is defined.

* #ifndef 'EXACTEXP'
*   if (count(acc,1)!=0);
*     argument accun;
*       id n=4-2*ep;
*     endargument;
*     id accun(x?)=acc(x);
*     repeat id acc(x1?)*acc(x2?) = acc(x1*x2);
*   endif;
*   #call ACCU(TabBN)
* ***if (count(accun,1)!=0) id acc(x?)=x;
* #endif

id agam=2;
id bgam=23/3;
id cgam=35/2+3*z2;
id dgam=275/192*16+23/2*z2-2*z3;
id egam=-16*(189/128 - 89/48*z3 - 3/32*z4 - 105/64*z2 - 9/64*z2^2);
id fgam = -384*(
         14917/18432 + 1/128*z3*z2 - 175/256*z3 + 649/1536*z4 + 1/320*z5
         - 275/3072*z2 - 23/1024*z2^2
         )
         +16*B4;

argument;
  id agam=2;
  id bgam=23/3;
  id cgam=35/2+3*z2;
  id dgam=275/192*16+23/2*z2-2*z3;
  id egam=-16*(189/128 - 89/48*z3 - 3/32*z4 - 105/64*z2 - 9/64*z2^2);
  id fgam = -384*(
         14917/18432 + 1/128*z3*z2 - 175/256*z3 + 649/1536*z4 + 1/320*z5
         - 275/3072*z2 - 23/1024*z2^2
         )
         +16*B4;
endargument;

#call ACCU()
        
* b int0,intt1;
* Print+s;
* .end

#endprocedure        


#procedure tad3l
*
* {top3l;6;3;0;1; ;(p1:3,1)(p2:4,3)(p3:3,2)(p4:1,4)(p5:4,2)(p6:2,1)}
* with all possible mass distributions, i.e.
* 
* {top3l;6;3;0;1; ;(p1:3,1)(p2:4,3)(p3:3,2)(p4:1,4)(p5:4,2)(p6:2,1);
*                  111111;011111;101111;110111;111011;111101;111110;
*                  001111;010111;011011;011101;011110;100111;101011;
*                  101101;101110;110011;110101;110110;111001;111010;
*                  111100;000111;001011;001101;001110;010011;010101;
*                  010110;011001;011010;011100;100011;100101;100110;
*                  101001;101010;101100;110001;110010;110100;111000;
*                  000011;000101;000110;001001;001010;001100;010001;
*                  010010;010100;011000;100001;100010;100100;101000;
*                  110000;000001;000010;000100;001000;010000;100000;
*                  000000\}

*
* do partial fractioning for all lines
*

        #do i = 1, 6
                #call partfrac(p`i',s`i'm)
        #enddo
        .sort

*
* discard massless tadpoles
*

        if (count(s1m,1,s2m,1,s3m,1,s4m,1,s5m,1,s6m,1)==0) discard;
        .sort

*
* map onto MATAD master topologies; mapping computed with EXP
* there are 2^6 - 1 = 63 cases
* (we dont care about symmetry ... )
*

        if ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,intd6;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p6, p2,p1, p3,p4, p4,-p3, p5,-p2, p6,p5,
        s2m,s1m, s3m,s4m, s4m,s3m, s5m,s2m, s6m,s5m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p6, p3,-p4, p4,-p3, p5,p5, p6,-p2,
        s1m,s1m, s3m,s4m, s4m,s3m, s5m,s5m, s6m,s2m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p6, p4,p2, p5,-p5, p6,p3, 
        s1m,s1m, s2m,s4m, s4m,s2m, s5m,s5m, s6m,s3m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p2, p4,-p6, p5,-p5, p6,-p4, 
        s1m,s1m, s2m,s3m, s3m,s2m, s5m,s5m, s6m,s4m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p2, p2,p1, p3,-p3, p4,p4, p5,p6, p6,p5, 
        s1m,s2m, s2m,s1m, s3m,s3m, s4m,s4m, s6m,s5m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p2, p3,p3, p4,p4, p5,p5, p6,p6, 
        s1m,s1m, s2m,s2m, s3m,s3m, s4m,s4m, s5m,s5m);
        multiply,intd5;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p3, p2,-p5, p3,-p6, p4,p2, p5,p4, p6,p1, 
        s3m,s6m, s4m,s2m, s5m,s4m, s6m,s1m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,-p6, p3,-p5, p4,p1, p5,p4, p6,p2, 
        s2m,s6m, s4m,s1m, s5m,s4m, s6m,s2m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p3, p2,p2, p3,p1, p4,-p5, p5,-p4, p6,-p6, 
        s2m,s2m, s3m,s1m, s5m,s4m, s6m,s6m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s2m,s3m, s3m,s6m, s4m,s5m, s6m,s4m);
        multiply,intbn;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p2, p4,-p6, p5,-p4, p6,-p5, 
        s2m,s1m, s3m,s2m, s4m,s6m, s5m,s4m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p6, p2,-p3, p3,p5, p4,p1, p5,p2, p6,p4, 
        s1m,s6m, s4m,s1m, s5m,s2m, s6m,s4m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p2, p2,p3, p3,-p1, p4,-p5, p5,-p6, p6,-p4, 
        s1m,s2m, s3m,s1m, s5m,s6m, s6m,s4m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p2, p4,-p6, p5,-p5, p6,-p4, 
        s1m,s1m, s3m,s2m, s4m,s6m, s6m,s4m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p6, p4,-p5, p5,-p4, p6,-p2, 
        s1m,s3m, s3m,s6m, s4m,s5m, s5m,s4m);
        multiply,intbn;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,-p5, p3,-p2, p4,p1, p5,p4, p6,p6, 
        s1m,s3m, s2m,s5m, s5m,s4m, s6m,s6m);
        multiply,intbn;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p2, p3,p3, p4,p4, p5,p5, p6,p6, 
        s1m,s1m, s2m,s2m, s4m,s4m, s6m,s6m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p2, p2,p1, p3,-p3, p4,p4, p5,p6, p6,p5, 
        s1m,s2m, s2m,s1m, s4m,s4m, s5m,s6m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p6, p3,-p4, p4,-p3, p5,p5, p6,-p2, 
        s1m,s1m, s2m,s6m, s3m,s4m, s6m,s2m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p6, p2,p1, p3,p4, p4,-p3, p5,-p2, p6,p5, 
        s1m,s6m, s2m,s1m, s3m,s4m, s5m,s2m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p6, p4,p2, p5,-p5, p6,p3, 
        s1m,s1m, s2m,s4m, s3m,s6m, s4m,s2m);
        multiply,intd4;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p2, p3,p3, p4,p4, p5,p5, p6,p6, 
        s4m,s4m, s5m,s5m, s6m,s6m);
        multiply,intbm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p4, p2,p5, p3,-p2, p4,p6, p5,p3, p6,-p1, 
        s3m,s2m, s5m,s3m, s6m,s1m);
        multiply,intdm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p5, p3,-p4, p4,-p3, p5,p2, p6,-p6, 
        s3m,s4m, s4m,s3m, s6m,s6m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p5, p2,p1, p3,p4, p4,-p3, p5,-p6, p6,p2, 
        s3m,s4m, s4m,s3m, s5m,s6m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p5, p2,-p3, p3,p2, p4,p1, p5,p6, p6,p4, 
        s2m,s3m, s5m,s6m, s6m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p5, p4,p6, p5,-p2, p6,p3, 
        s2m,s4m, s4m,s6m, s6m,s3m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p1, p3,-p6, p4,p2, p5,p3, p6,-p5, 
        s2m,s1m, s4m,s2m, s5m,s3m);
        multiply,intdm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s2m,s3m, s3m,s6m, s6m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p6, p4,p2, p5,-p5, p6,p3, 
        s2m,s4m, s3m,s6m, s5m,s5m);
        multiply,intbm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p6, p3,p3, p4,p4, p5,p2, p6,p5, 
        s2m,s6m, s3m,s3m, s4m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,-p5, p3,-p2, p4,p1, p5,p4, p6,p6, 
        s1m,s3m, s5m,s4m, s6m,s6m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p6, p4,p2, p5,-p5, p6,p3,  
        s1m,s1m, s4m,s2m, s6m,s3m);
        multiply,intdm;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p1, p3,-p5, p4,p6, p5,p3, p6,-p2, 
        s1m,s4m, s4m,s6m, s5m,s3m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p4, p2,p1, p3,-p6, p4,p2, p5,p3, p6,-p5, 
        s1m,s4m, s3m,s6m, s6m,s5m);
        multiply,intbm;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p6, p4,-p5, p5,-p4, p6,-p2, 
        s1m,s3m, s3m,s6m, s5m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p6, p2,p1, p3,-p3, p4,p4, p5,p5, p6,p2, 
        s1m,s6m, s3m,s3m, s4m,s4m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p6, p2,p4, p3,p2, p4,p1, p5,-p5, p6,-p3, 
        s1m,s6m, s2m,s4m, s6m,s3m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p6, p3,-p2, p4,p1, p5,-p3, p6,-p5, 
        s1m,s4m, s2m,s6m, s5m,s3m);
        multiply,intbn1;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p5, p3,-p2, p4,p6, p5,p3, p6,-p1, 
        s1m,s4m, s2m,s5m, s4m,s6m);
        multiply,intbm;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p2, p3,p3, p4,p4, p5,p5, p6,p6, 
        s1m,s1m, s2m,s2m, s3m,s3m);
        multiply,intdm;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,-p5, p3,-p2, p4,p1, p5,p4, p6,p6, 
        s5m,s4m, s6m,s6m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p2, p2,-p3, p3,-p5, p4,-p6, p5,-p1, p6,-p4, 
        s4m,s6m, s6m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p2, p3,p5, p4,-p6, p5,-p4, p6,-p1, 
        s4m,s6m, s5m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s3m,s6m, s6m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p6, p4,-p5, p5,-p4, p6,-p2, 
        s3m,s6m, s5m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,-p5, p3,-p2, p4,p1, p5,p4, p6,p6, 
        s3m,s2m, s4m,s1m);
        multiply,intdn;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p3, p2,p1, p3,p6, p4,-p5, p5,-p4, p6,-p2, 
        s2m,s1m, s6m,s2m);
        multiply,intdn;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p3, p2,-p6, p3,-p1, p4,p2, p5,p4, p6,p5, 
        s2m,s6m, s5m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p5, p4,p6, p5,-p2, p6,p3, 
        s2m,s4m, s4m,s6m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p2, p2,-p6, p3,-p4, p4,-p3, p5,p1, p6,-p5, 
        s2m,s6m, s3m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,-p6, p2,-p3, p3,p1, p4,p2, p5,p5, p6,p4, 
        s1m,s6m, s6m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s1m,s1m, s5m,s2m);
        multiply,intdn;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p1, p3,-p5, p4,p6, p5,p3, p6,-p2, 
        s1m,s4m, s4m,s6m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p6, p2,p2, p3,p4, p4,-p3, p5,-p5, p6,p1, 
        s1m,s6m, s3m,s4m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p4, p2,p6, p3,-p2, p4,p1, p5,-p3, p6,-p5, 
        s1m,s4m, s2m,s6m);
        multiply,intbn2;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)!=0) );
*
        multiply,replace_(p1,p1, p2,-p5, p3,-p4, p4,-p3, p5,p2, p6,-p6, 
        s6m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)!=0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,-p5, p2,p1, p3,p4, p4,-p3, p5,-p6, p6,p2, 
        s5m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)!=0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p4, p3,p5, p4,p6, p5,-p2, p6,p3, 
        s4m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)==0) && (count(s3m,1)!=0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,-p3, p3,-p6, p4,-p5, p5,-p2, p6,-p4, 
        s3m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)==0) && (count(s2m,1)!=0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p1, p2,p6, p3,p3, p4,p4, p5,p2, p6,p5, 
        s2m,s6m);
        multiply,intbn3;
*
        elseif ( (count(s1m,1)!=0) && (count(s2m,1)==0) && (count(s3m,1)==0) 
        && (count(s4m,1)==0) && (count(s5m,1)==0) && (count(s6m,1)==0) );
*
        multiply,replace_(p1,p6, p2,p1, p3,-p3, p4,p4, p5,p5, p6,p2, 
        s1m,s6m);
        multiply,intbn3;
        endif;
        .sort

*
* apply recurrence relations
* 

        #call dorec3l
        

        #call Conv2exact
        #call DoG
        #call subSimple
        #call GammaArgToOne
        
#endprocedure



#procedure matad(LOOPS)
        #message MATAD called for `LOOPS'-loop integral

        #call tad`LOOPS'l
        

*         b intn1,x6;        
*         Print+s;        
*         .end

*         if(count(int0,1));
*         Multiply 1/int0;  
        
*         else;
*         exit "Not all integrals reduced";

*         endif;        
#endprocedure


#procedure subSimple

* 
* Here we substitute all integrals known in terms of Gamma functions
* 
*   gm3norm(a1,a2,a3) - Two loop tadpole    [MM0]
*   gm2norm(a1,a2)    - One-loop tadpole    [M0]
*   G(a1,a2)          - One-loop propagator [00]
*         


*       Upto ep^5        
*           (M,M,0)   
*         id gm3norm(0,0,0) = -7 - ep^(-2) - 3/ep - z2 + ep*(-15 - 3*z2 + (2*z3)/3) +
*         ep^2*(-31 - 7*z2 + 2*z3 - (7*z4)/4) +
*         ep^3*(-63 - 15*z2 + 12*z2^2 + (14*z3)/3 + (2*z2*z3)/3 - (141*z4)/4 + (2*z5)/5) +
*         ep^4*(-127 - 31*z2 + 28*z2^2 + (17*z2^3)/4 + 10*z3 + 2*z2*z3 - (2*z3^2)/9 -
*         (329*z4)/4 + (81*z2*z4)/16 + (6*z5)/5 - (1881*z6)/64) +
*         ep^5*(-255 - 63*z2 + 60*z2^2 + (51*z2^3)/4 + (62*z3)/3 + (14*z2*z3)/3 -
*         (8*z2^2*z3)/3 - (2*z3^2)/3 - (705*z4)/4 + (243*z2*z4)/16 + (47*z3*z4)/6 +
*         (14*z5)/5 + (2*z2*z5)/5 - (5643*z6)/64 + (2*z7)/7);
        
*         id gm3norm(0,0,1) = -13/3 - 1/(3*ep^2) - 4/(3*ep) - (2*z2)/3 + ep*(-40/3 - (8*z2)/3 - (28*z3)/9) +
*         ep^2*(-121/3 - (26*z2)/3 - (112*z3)/9 + 9*z4) +
*         ep^3*(-364/3 - (80*z2)/3 + (185*z2^2)/3 - (364*z3)/9 - (56*z2*z3)/9 - (709*z4)/6 -
*         (748*z5)/15) + ep^4*(-1093/3 - (242*z2)/3 + (2405*z2^2)/12 + (455*z2^3)/6 -
*         (1120*z3)/9 - (224*z2*z3)/9 - (392*z3^2)/27 - (9217*z4)/24 - (185*z2*z4)/12 -
*         (2992*z5)/15 - (2231*z6)/24) + ep^5*(-3280/3 - (728*z2)/3 + (1850*z2^2)/3 +
*         (910*z2^3)/3 - (3388*z3)/9 - (728*z2*z3)/9 + (1295*z2^2*z3)/9 - (1568*z3^2)/27 -
*         (3545*z4)/3 - (185*z2*z4)/3 - (4963*z3*z4)/18 - (9724*z5)/15 - (1496*z2*z5)/15 -
*         (2231*z6)/6 - (14068*z7)/21);

*         id gm3norm(0,0,2) = -7/2 - 1/(6*ep^2) - 5/(6*ep) - z2/2 + ep*(-85/6 - (5*z2)/2 - (62*z3)/9) +
*         ep^2*(-341/6 - (21*z2)/2 - (310*z3)/9 + (251*z4)/8) +
*         ep^3*(-455/2 - (85*z2)/2 + (550*z2^2)/3 - (434*z3)/3 - (62*z2*z3)/3 - (7235*z4)/24 -
*         (3254*z5)/15) + ep^4*(-5461/6 - (341*z2)/2 + 770*z2^2 + (3135*z2^3)/8 -
*         (5270*z3)/9 - (310*z2*z3)/3 - (3844*z3^2)/27 - (10129*z4)/8 - (12419*z2*z4)/96 -
*         (3254*z5)/3 - (74989*z6)/384) + ep^5*(-21845/6 - (1365*z2)/2 + (9350*z2^2)/3 +
*         (15675*z2^3)/8 - (21142*z3)/9 - 434*z2*z3 + (13640*z2^2*z3)/9 - (19220*z3^2)/27 -
*         (122995*z4)/24 - (62095*z2*z4)/96 - (44857*z3*z4)/18 - (22778*z5)/5 -
*         (3254*z2*z5)/5 - (374945*z6)/384 - (130682*z7)/21);

*         id gm3norm(0,0,3) = -31/10 - 1/(10*ep^2) - 3/(5*ep) - (2*z2)/5 + ep*(-78/5 - (12*z2)/5 - (161*z3)/15) +
*         ep^2*(-781/10 - (62*z2)/5 - (322*z3)/5 + (328*z4)/5) +
*         ep^3*(-1953/5 - (312*z2)/5 + (8139*z2^2)/20 - (4991*z3)/15 - (644*z2*z3)/15 -
*         (24951*z4)/40 - (14309*z5)/25) + ep^4*(-19531/10 - (1562*z2)/5 + (84103*z2^2)/40 +
*         (6247*z2^3)/5 - (8372*z3)/5 - (1288*z2*z3)/5 - (25921*z3^2)/45 - (257827*z4)/80 -
*         450*z2*z4 - (85854*z5)/25 - (12249*z6)/40) +
*         ep^5*(-48828/5 - (7812*z2)/5 + (105807*z2^2)/10 + (37482*z2^3)/5 - (125741*z3)/15 -
*         (19964*z2*z3)/15 + (436793*z2^2*z3)/60 - (51842*z3^2)/15 - (324363*z4)/20 -
*         2700*z2*z4 - (1339037*z3*z4)/120 - (443579*z5)/25 - (57236*z2*z5)/25 -
*         (36747*z6)/20 - (1001321*z7)/35);

*         id gm3norm(0,0,4) = -43/15 - 1/(15*ep^2) - 7/(15*ep) - z2/3 + ep*(-259/15 - (7*z2)/3 - (658*z3)/45) +
*         ep^2*(-311/3 - (43*z2)/3 - (4606*z3)/45 + (447*z4)/4) +
*         ep^3*(-9331/15 - (259*z2)/3 + (11438*z2^2)/15 - (28294*z3)/45 - (658*z2*z3)/9 -
*         (13489*z4)/12 - (89098*z5)/75) + ep^4*(-55987/15 - (1555*z2)/3 + (70262*z2^2)/15 +
*         (36785*z2^3)/12 - (170422*z3)/45 - (4606*z2*z3)/9 - (216482*z3^2)/135 -
*         (82861*z4)/12 - (271811*z2*z4)/240 - (623686*z5)/75 - (352693*z6)/960) +
*         ep^5*(-335923/15 - (9331*z2)/3 + (423206*z2^2)/15 + (257495*z2^3)/12 - (204638*z3)/9 -
*         (28294*z2*z3)/9 + (1075172*z2^2*z3)/45 - (1515374*z3^2)/135 - (499093*z4)/12 -
*         (1902677*z2*z4)/240 - (633983*z3*z4)/18 - (3831214*z5)/75 - (89098*z2*z5)/15 -
*         (2468851*z6)/960 - (9641938*z7)/105);

* *       Upto ep^10        
*         id gm2norm(0,0) = 1 + ep + ep^2*(1 + z2/2) + ep^3*(1 + z2/2 - z3/3) +
*         ep^4*(1 + z2/2 - z2^2/4 - z3/3 + (19*z4)/16) +
*         ep^5*(1 + z2/2 - z2^2/4 - z3/3 - (z2*z3)/6 + (19*z4)/16 - z5/5) +
*         ep^6*(1 + z2/2 - z2^2/4 - z2^3/8 - z3/3 - (z2*z3)/6 + z3^2/18 + (19*z4)/16 +
*         (z2*z4)/16 - z5/5 + (117*z6)/128) +
*         ep^7*(1 + z2/2 - z2^2/4 - z2^3/8 - z3/3 - (z2*z3)/6 + (z2^2*z3)/12 + z3^2/18 +
*         (19*z4)/16 + (z2*z4)/16 - (19*z3*z4)/48 - z5/5 - (z2*z5)/10 + (117*z6)/128 -
*         z7/7) + ep^8*(1 + z2/2 - z2^2/4 - z2^3/8 - z3/3 - (z2*z3)/6 + (z2^2*z3)/12 +
*         z3^2/18 + (z2*z3^2)/36 + (19*z4)/16 + (z2*z4)/16 - (17*z2^2*z4)/64 -
*         (19*z3*z4)/48 + (99*z4^2)/256 - z5/5 - (z2*z5)/10 + (z3*z5)/15 +
*         (117*z6)/128 - (5*z2*z6)/128 - z7/7 + (2455*z8)/3072) +
*         ep^9*(1 + z2/2 - z2^2/4 - z2^3/8 - z3/3 - (z2*z3)/6 + (z2^2*z3)/12 +
*         (z2^3*z3)/24 + z3^2/18 + (z2*z3^2)/36 - z3^3/162 + (19*z4)/16 + (z2*z4)/16 -
*         (17*z2^2*z4)/64 - (19*z3*z4)/48 - (z2*z3*z4)/48 + (99*z4^2)/256 - z5/5 -
*         (z2*z5)/10 + (z2^2*z5)/20 + (z3*z5)/15 - (19*z4*z5)/80 + (117*z6)/128 -
*         (5*z2*z6)/128 - (39*z3*z6)/128 - z7/7 - (z2*z7)/14 + (2455*z8)/3072 -
*         z9/9) + ep^10*(1 + (12883*z10)/20480 + z2/2 - z2^2/4 - z2^3/8 + z2^5/960 -
*         z3/3 - (z2*z3)/6 + (z2^2*z3)/12 + (z2^3*z3)/24 + z3^2/18 + (z2*z3^2)/36 -
*         (z2^2*z3^2)/72 - z3^3/162 + (19*z4)/16 + (z2*z4)/16 - (17*z2^2*z4)/64 +
*         (7*z2^3*z4)/768 - (19*z3*z4)/48 - (z2*z3*z4)/48 + (19*z3^2*z4)/288 +
*         (99*z4^2)/256 - (67*z2*z4^2)/512 - z5/5 - (z2*z5)/10 + (z2^2*z5)/20 +
*         (z3*z5)/15 + (z2*z3*z5)/30 - (19*z4*z5)/80 + z5^2/50 + (117*z6)/128 -
*         (5*z2*z6)/128 - (683*z2^2*z6)/3072 - (39*z3*z6)/128 + (1565*z4*z6)/3072 -
*         z7/7 - (z2*z7)/14 + (z3*z7)/21 + (2455*z8)/3072 - (67*z2*z8)/6144 - z9/9);

* * Insertions into the massless line
*         id gm2norm(0,1) = 1 + ep + ep^2*(1 + (7*z2)/2) + ep^3*(1 + (7*z2)/2 - z3/3) +
*         ep^4*(1 + (7*z2)/2 + 2*z2^2 - z3/3 + (289*z4)/16) +
*         ep^5*(1 + (7*z2)/2 + 2*z2^2 - z3/3 - (7*z2*z3)/6 + (289*z4)/16 - z5/5) +
*         ep^6*(1 + (7*z2)/2 + 2*z2^2 - 2*z2^3 - z3/3 - (7*z2*z3)/6 + z3^2/18 + (289*z4)/16 +
*         (109*z2*z4)/4 - z5/5 + (7803*z6)/128) +
*         ep^7*(1 + (7*z2)/2 + 2*z2^2 - 2*z2^3 - z3/3 - (7*z2*z3)/6 - (2*z2^2*z3)/3 + z3^2/18 +
*         (289*z4)/16 + (109*z2*z4)/4 - (289*z3*z4)/48 - z5/5 - (7*z2*z5)/10 + (7803*z6)/128 -
*         z7/7) + ep^8*(1 + (7*z2)/2 + 2*z2^2 - 2*z2^3 - z3/3 - (7*z2*z3)/6 - (2*z2^2*z3)/3 +
*         z3^2/18 + (7*z2*z3^2)/36 + (289*z4)/16 + (109*z2*z4)/4 - (71*z2^2*z4)/4 -
*         (289*z3*z4)/48 + (657*z4^2)/8 - z5/5 - (7*z2*z5)/10 + (z3*z5)/15 + (7803*z6)/128 +
*         (2923*z2*z6)/32 - z7/7 + (645565*z8)/3072) +
*         ep^9*(1 + (7*z2)/2 + 2*z2^2 - 2*z2^3 - z3/3 - (7*z2*z3)/6 - (2*z2^2*z3)/3 +
*         (2*z2^3*z3)/3 + z3^2/18 + (7*z2*z3^2)/36 - z3^3/162 + (289*z4)/16 + (109*z2*z4)/4 -
*         (71*z2^2*z4)/4 - (289*z3*z4)/48 - (109*z2*z3*z4)/12 + (657*z4^2)/8 - z5/5 -
*         (7*z2*z5)/10 - (2*z2^2*z5)/5 + (z3*z5)/15 - (289*z4*z5)/80 + (7803*z6)/128 +
*         (2923*z2*z6)/32 - (2601*z3*z6)/128 - z7/7 - (z2*z7)/2 + (645565*z8)/3072 - z9/9) +
*         ep^10*(1 + (14929*z10)/20 + (7*z2)/2 + 2*z2^2 - 2*z2^3 + z2^5/960 - z3/3 -
*         (7*z2*z3)/6 - (2*z2^2*z3)/3 + (2*z2^3*z3)/3 + z3^2/18 + (7*z2*z3^2)/36 +
*         (z2^2*z3^2)/9 - z3^3/162 + (289*z4)/16 + (109*z2*z4)/4 - (71*z2^2*z4)/4 +
*         (7*z2^3*z4)/768 - (289*z3*z4)/48 - (109*z2*z3*z4)/12 + (289*z3^2*z4)/288 +
*         (657*z4^2)/8 - (4897*z2*z4^2)/128 - z5/5 - (7*z2*z5)/10 - (2*z2^2*z5)/5 +
*         (z3*z5)/15 + (7*z2*z3*z5)/30 - (289*z4*z5)/80 + z5^2/50 + (7803*z6)/128 +
*         (2923*z2*z6)/32 - (187793*z2^2*z6)/3072 - (2601*z3*z6)/128 + (847975*z4*z6)/1536 -
*         z7/7 - (z2*z7)/2 + (z3*z7)/21 + (645565*z8)/3072 + (1936427*z2*z8)/6144 - z9/9);

*         id gm2norm(0,2) =1 + ep + ep^2*(1 + (17*z2)/2) + ep^3*(1 + (17*z2)/2 - z3/3) +
*         ep^4*(1 + (17*z2)/2 + (63*z2^2)/4 - z3/3 + (1459*z4)/16) +
*         ep^5*(1 + (17*z2)/2 + (63*z2^2)/4 - z3/3 - (17*z2*z3)/6 + (1459*z4)/16 - z5/5) +
*         ep^6*(1 + (17*z2)/2 + (63*z2^2)/4 - (81*z2^3)/8 - z3/3 - (17*z2*z3)/6 + z3^2/18 +
*         (1459*z4)/16 + (5841*z2*z4)/16 - z5/5 + (88933*z6)/128) +
*         ep^7*(1 + (17*z2)/2 + (63*z2^2)/4 - (81*z2^3)/8 - z3/3 - (17*z2*z3)/6 -
*         (21*z2^2*z3)/4 + z3^2/18 + (1459*z4)/16 + (5841*z2*z4)/16 - (1459*z3*z4)/48 - z5/5 -
*         (17*z2*z5)/10 + (88933*z6)/128 - z7/7) +
*         ep^8*(1 + (17*z2)/2 + (63*z2^2)/4 - (81*z2^3)/8 - z3/3 - (17*z2*z3)/6 -
*         (21*z2^2*z3)/4 + z3^2/18 + (17*z2*z3^2)/36 + (1459*z4)/16 + (5841*z2*z4)/16 -
*         (13041*z2^2*z4)/64 - (1459*z3*z4)/48 + (532899*z4^2)/256 - z5/5 - (17*z2*z5)/10 +
*         (z3*z5)/15 + (88933*z6)/128 + (355707*z2*z6)/128 - z7/7 + (16546775*z8)/3072) +
*         ep^9*(1 + (17*z2)/2 + (63*z2^2)/4 - (81*z2^3)/8 - z3/3 - (17*z2*z3)/6 -
*         (21*z2^2*z3)/4 + (27*z2^3*z3)/8 + z3^2/18 + (17*z2*z3^2)/36 - z3^3/162 +
*         (1459*z4)/16 + (5841*z2*z4)/16 - (13041*z2^2*z4)/64 - (1459*z3*z4)/48 -
*         (1947*z2*z3*z4)/16 + (532899*z4^2)/256 - z5/5 - (17*z2*z5)/10 - (63*z2^2*z5)/20 +
*         (z3*z5)/15 - (1459*z4*z5)/80 + (88933*z6)/128 + (355707*z2*z6)/128 -
*         (88933*z3*z6)/384 - z7/7 - (17*z2*z7)/14 + (16546775*z8)/3072 - z9/9) +
*         ep^10*(1 + (881658571*z10)/20480 + (17*z2)/2 + (63*z2^2)/4 - (81*z2^3)/8 + z2^5/960 -
*         z3/3 - (17*z2*z3)/6 - (21*z2^2*z3)/4 + (27*z2^3*z3)/8 + z3^2/18 + (17*z2*z3^2)/36 +
*         (7*z2^2*z3^2)/8 - z3^3/162 + (1459*z4)/16 + (5841*z2*z4)/16 - (13041*z2^2*z4)/64 +
*         (7*z2^3*z4)/768 - (1459*z3*z4)/48 - (1947*z2*z3*z4)/16 + (1459*z3^2*z4)/288 +
*         (532899*z4^2)/256 - (518323*z2*z4^2)/512 - z5/5 - (17*z2*z5)/10 - (63*z2^2*z5)/20 +
*         (z3*z5)/15 + (17*z2*z3*z5)/30 - (1459*z4*z5)/80 + z5^2/50 + (88933*z6)/128 +
*         (355707*z2*z6)/128 - (4805003*z2^2*z6)/3072 - (88933*z3*z6)/384 +
*         (97376045*z4*z6)/3072 - z7/7 - (17*z2*z7)/14 + (z3*z7)/21 + (16546775*z8)/3072 +
*        (132373597*z2*z8)/6144 - z9/9);

*         id gm2norm(0,3) =1 + ep + ep^2*(1 + (31*z2)/2) + ep^3*(1 + (31*z2)/2 - z3/3) +
*         ep^4*(1 + (31*z2)/2 + 56*z2^2 - z3/3 + (4609*z4)/16) +
*         ep^5*(1 + (31*z2)/2 + 56*z2^2 - z3/3 - (31*z2*z3)/6 + (4609*z4)/16 - z5/5) +
*         ep^6*(1 + (31*z2)/2 + 56*z2^2 - 32*z2^3 - z3/3 - (31*z2*z3)/6 + z3^2/18 +
*         (4609*z4)/16 + 2161*z2*z4 - z5/5 + (499707*z6)/128) +
*         ep^7*(1 + (31*z2)/2 + 56*z2^2 - 32*z2^3 - z3/3 - (31*z2*z3)/6 - (56*z2^2*z3)/3 +
*         z3^2/18 + (4609*z4)/16 + 2161*z2*z4 - (4609*z3*z4)/48 - z5/5 - (31*z2*z5)/10 +
*         (499707*z6)/128 - z7/7) + ep^8*(1 + (31*z2)/2 + 56*z2^2 - 32*z2^3 - z3/3 -
*         (31*z2*z3)/6 - (56*z2^2*z3)/3 + z3^2/18 + (31*z2*z3^2)/36 + (4609*z4)/16 +
*         2161*z2*z4 - 1148*z2^2*z4 - (4609*z3*z4)/48 + 20754*z4^2 - z5/5 - (31*z2*z5)/10 +
*         (z3*z5)/15 + (499707*z6)/128 + (234235*z2*z6)/8 - z7/7 + (165281725*z8)/3072) +
*         ep^9*(1 + (31*z2)/2 + 56*z2^2 - 32*z2^3 - z3/3 - (31*z2*z3)/6 - (56*z2^2*z3)/3 +
*         (32*z2^3*z3)/3 + z3^2/18 + (31*z2*z3^2)/36 - z3^3/162 + (4609*z4)/16 + 2161*z2*z4 -
*         1148*z2^2*z4 - (4609*z3*z4)/48 - (2161*z2*z3*z4)/3 + 20754*z4^2 - z5/5 -
*         (31*z2*z5)/10 - (56*z2^2*z5)/5 + (z3*z5)/15 - (4609*z4*z5)/80 + (499707*z6)/128 +
*         (234235*z2*z6)/8 - (166569*z3*z6)/128 - z7/7 - (31*z2*z7)/14 + (165281725*z8)/3072 -
*         z9/9) + ep^10*(1 + (7644671*z10)/10 + (31*z2)/2 + 56*z2^2 - 32*z2^3 + z2^5/960 -
*         z3/3 - (31*z2*z3)/6 - (56*z2^2*z3)/3 + (32*z2^3*z3)/3 + z3^2/18 + (31*z2*z3^2)/36 +
*         (28*z2^2*z3^2)/9 - z3^3/162 + (4609*z4)/16 + 2161*z2*z4 - 1148*z2^2*z4 +
*         (7*z2^3*z4)/768 - (4609*z3*z4)/48 - (2161*z2*z3*z4)/3 + (4609*z3^2*z4)/288 +
*         20754*z4^2 - (1308673*z2*z4^2)/128 - z5/5 - (31*z2*z5)/10 - (56*z2^2*z5)/5 +
*         (z3*z5)/15 + (31*z2*z3*z5)/30 - (4609*z4*z5)/80 + z5^2/50 + (499707*z6)/128 +
*         (234235*z2*z6)/8 - (47979953*z2^2*z6)/3072 - (166569*z3*z6)/128 +
*         (863859775*z4*z6)/1536 - z7/7 - (31*z2*z7)/14 + (z3*z7)/21 + (165281725*z8)/3072 +
*         (2479224803*z2*z8)/6144 - z9/9);

*         id gm2norm(0,4) =1 + ep + ep^2*(1 + (49*z2)/2) + ep^3*(1 + (49*z2)/2 - z3/3) +
*         ep^4*(1 + (49*z2)/2 + (575*z2^2)/4 - z3/3 + (11251*z4)/16) +
*         ep^5*(1 + (49*z2)/2 + (575*z2^2)/4 - z3/3 - (49*z2*z3)/6 + (11251*z4)/16 - z5/5) +
*         ep^6*(1 + (49*z2)/2 + (575*z2^2)/4 - (625*z2^3)/8 - z3/3 - (49*z2*z3)/6 + z3^2/18 +
*         (11251*z4)/16 + (135025*z2*z4)/16 - z5/5 + (1906245*z6)/128) +
*         ep^7*(1 + (49*z2)/2 + (575*z2^2)/4 - (625*z2^3)/8 - z3/3 - (49*z2*z3)/6 -
*         (575*z2^2*z3)/12 + z3^2/18 + (11251*z4)/16 + (135025*z2*z4)/16 - (11251*z3*z4)/48 -
*         z5/5 - (49*z2*z5)/10 + (1906245*z6)/128 - z7/7) +
*         ep^8*(1 + (49*z2)/2 + (575*z2^2)/4 - (625*z2^3)/8 - z3/3 - (49*z2*z3)/6 -
*         (575*z2^2*z3)/12 + z3^2/18 + (49*z2*z3^2)/36 + (11251*z4)/16 + (135025*z2*z4)/16 -
*         (280625*z2^2*z4)/64 - (11251*z3*z4)/48 + (31651875*z4^2)/256 - z5/5 - (49*z2*z5)/10 +
*         (z3*z5)/15 + (1906245*z6)/128 + (22874875*z2*z6)/128 - z7/7 + (985156183*z8)/3072) +
*         ep^9*(1 + (49*z2)/2 + (575*z2^2)/4 - (625*z2^3)/8 - z3/3 - (49*z2*z3)/6 -
*         (575*z2^2*z3)/12 + (625*z2^3*z3)/24 + z3^2/18 + (49*z2*z3^2)/36 - z3^3/162 +
*         (11251*z4)/16 + (135025*z2*z4)/16 - (280625*z2^2*z4)/64 - (11251*z3*z4)/48 -
*         (135025*z2*z3*z4)/48 + (31651875*z4^2)/256 - z5/5 - (49*z2*z5)/10 - (115*z2^2*z5)/4 +
*         (z3*z5)/15 - (11251*z4*z5)/80 + (1906245*z6)/128 + (22874875*z2*z6)/128 -
*         (635415*z3*z6)/128 - z7/7 - (7*z2*z7)/2 + (985156183*z8)/3072 - z9/9) +
*         ep^10*(1 + (145810544827*z10)/20480 + (49*z2)/2 + (575*z2^2)/4 - (625*z2^3)/8 +
*         z2^5/960 - z3/3 - (49*z2*z3)/6 - (575*z2^2*z3)/12 + (625*z2^3*z3)/24 + z3^2/18 +
*         (49*z2*z3^2)/36 + (575*z2^2*z3^2)/72 - z3^3/162 + (11251*z4)/16 + (135025*z2*z4)/16 -
*         (280625*z2^2*z4)/64 + (7*z2^3*z4)/768 - (11251*z3*z4)/48 - (135025*z2*z3*z4)/48 +
*         (11251*z3^2*z4)/288 + (31651875*z4^2)/256 - (31359379*z2*z4^2)/512 - z5/5 -
*         (49*z2*z5)/10 - (115*z2^2*z5)/4 + (z3*z5)/15 + (49*z2*z3*z5)/30 - (11251*z4*z5)/80 +
*         z5^2/50 + (1906245*z6)/128 + (22874875*z2*z6)/128 - (285956171*z2^2*z6)/3072 -
*         (635415*z3*z6)/128 + (16086759245*z4*z6)/3072 - z7/7 - (7*z2*z7)/2 + (z3*z7)/21 +
*         (985156183*z8)/3072 + (23643746717*z2*z8)/6144 - z9/9);

*         id gm2norm(0,5) =1 + ep + ep^2*(1 + (71*z2)/2) + ep^3*(1 + (71*z2)/2 - z3/3) +
*         ep^4*(1 + (71*z2)/2 + 306*z2^2 - z3/3 + (23329*z4)/16) +
*         ep^5*(1 + (71*z2)/2 + 306*z2^2 - z3/3 - (71*z2*z3)/6 + (23329*z4)/16 - z5/5) +
*         ep^6*(1 + (71*z2)/2 + 306*z2^2 - 162*z2^3 - z3/3 - (71*z2*z3)/6 + z3^2/18 +
*         (23329*z4)/16 + (102069*z2*z4)/4 - z5/5 + (5692027*z6)/128) +
*         ep^7*(1 + (71*z2)/2 + 306*z2^2 - 162*z2^3 - z3/3 - (71*z2*z3)/6 - 102*z2^2*z3 +
*         z3^2/18 + (23329*z4)/16 + (102069*z2*z4)/4 - (23329*z3*z4)/48 - z5/5 -
*         (71*z2*z5)/10 + (5692027*z6)/128 - z7/7) +
*         ep^8*(1 + (71*z2)/2 + 306*z2^2 - 162*z2^3 - z3/3 - (71*z2*z3)/6 - 102*z2^2*z3 +
*         z3^2/18 + (71*z2*z3^2)/36 + (23329*z4)/16 + (102069*z2*z4)/4 - (52407*z2^2*z4)/4 -
*         (23329*z3*z4)/48 + (4252257*z4^2)/8 - z5/5 - (71*z2*z5)/10 + (z3*z5)/15 +
*         (5692027*z6)/128 + (24902595*z2*z6)/32 - z7/7 + (4235991485*z8)/3072) +
*         ep^9*(1 + (71*z2)/2 + 306*z2^2 - 162*z2^3 - z3/3 - (71*z2*z3)/6 - 102*z2^2*z3 +
*         54*z2^3*z3 + z3^2/18 + (71*z2*z3^2)/36 - z3^3/162 + (23329*z4)/16 +
*         (102069*z2*z4)/4 - (52407*z2^2*z4)/4 - (23329*z3*z4)/48 - (34023*z2*z3*z4)/4 +
*         (4252257*z4^2)/8 - z5/5 - (71*z2*z5)/10 - (306*z2^2*z5)/5 + (z3*z5)/15 -
*         (23329*z4*z5)/80 + (5692027*z6)/128 + (24902595*z2*z6)/32 - (5692027*z3*z6)/384 -
*         z7/7 - (71*z2*z7)/14 + (4235991485*z8)/3072 - z9/9) +
*         ep^10*(1 + (881660617*z10)/20 + (71*z2)/2 + 306*z2^2 - 162*z2^3 + z2^5/960 - z3/3 -
*         (71*z2*z3)/6 - 102*z2^2*z3 + 54*z2^3*z3 + z3^2/18 + (71*z2*z3^2)/36 + 17*z2^2*z3^2 -
*         z3^3/162 + (23329*z4)/16 + (102069*z2*z4)/4 - (52407*z2^2*z4)/4 + (7*z2^3*z4)/768 -
*         (23329*z3*z4)/48 - (34023*z2*z3*z4)/4 + (23329*z3^2*z4)/288 + (4252257*z4^2)/8 -
*         (33802273*z2*z4^2)/128 - z5/5 - (71*z2*z5)/10 - (306*z2^2*z5)/5 + (z3*z5)/15 +
*         (71*z2*z3*z5)/30 - (23329*z4*z5)/80 + z5^2/50 + (5692027*z6)/128 +
*         (24902595*z2*z6)/32 - (1229517713*z2^2*z6)/3072 - (5692027*z3*z6)/384 +
*         (49798077415*z4*z6)/1536 - z7/7 - (71*z2*z7)/14 + (z3*z7)/21 + (4235991485*z8)/3072 +
*         (148259699563*z2*z8)/6144 - z9/9);

* Insertions into the massive line        
*         id gm2norm(1,0) =;
*         id gm2norm(2,0) =;
*         id gm2norm(3,0) =;
*         id gm2norm(4,0) =;
*         id gm2norm(5,0) =;                


        id gm3norm(x1?,x2?,x3?) = 
        Gam(-1,2+x1+x2+x3)*
        Gam(0,1+x1+x3)*
        Gam(0,1+x2+x3)*
        Gam(1,-1-x3)*
        iGam(1,x1 )*
        iGam(1,x2 )*
        iGam(0,2+x1+x2+2*x3)*
        iGam(2,-1);        
        
        id gm2norm(k2?,k4?) =
        Gam(1,-1-k4)*
        Gam(1,1+k2+k4)*
        iGam(2,-1)*
        iGam(1,k2);


        id GschemeConstants(x1?,x2?)= 
        Gam(1,1+x1+x2)*
        Gam(1,-1-x1)*
        Gam(1,-1-x2)*
        iGam(1,x1)*
        iGam(1,x2)*
        iGam(2,-2-x1-x2)*
        rat(1,ep);
*         rat(1,ep^(1+x1+x2));
*       We normalize G function to have 1+a1+a2 maximum pole       



* id	GschemeConstants(2,0) = GschemeConstants(0,0)*(
* 		1
*       +ep*(2)
*       +ep^2*(8)
*       +ep^3*(32-16*z3)
*       +ep^4*(128-24*z4-32*z3)
*       +ep^5*(512-192*z5-48*z4-128*z3)
*       +ep^6*(2048-440*z6-384*z5-192*z4-512*z3+128*z3^2)
*       +ep^7*(8192-2304*z7-880*z6-1536*z5-768*z4-2048*z3+384*z3*z4+256*z3^2
* 		));
* id	GschemeConstants(1,0) = GschemeConstants(0,0)*(
* 		1
*       +ep*(1)
*       +ep^2*(3)
*       +ep^3*(9-6*z3)
*       +ep^4*(27-9*z4-6*z3)
*       +ep^5*(81-42*z5-9*z4-18*z3)
*       +ep^6*(243-90*z6-42*z5-27*z4-54*z3+18*z3^2)
*       +ep^7*(729-294*z7-90*z6-126*z5-81*z4-162*z3+54*z3*z4+18*z3^2
* 		));

* id	GschemeConstants(0,0)^3 = epp^3*(1
*       +ep*(6)
*       +ep^2*(24)
*       +ep^3*(80-7*z3)
*       +ep^4*(240-39/4*z4-42*z3)
*       +ep^5*(672-93/5*z5-117/2*z4-168*z3)
*       +ep^6*(1792-61/2*z6-558/5*z5-234*z4-560*z3+49/2*z3^2)
*       +ep^7*(4608-381/7*z7-183*z6-2232/5*z5-780*z4-1680*z3+273/4*z3*z4+
*          147*z3^2));
* id	GschemeConstants(0,0)^2 = epp^2*(1
*       +ep*(4)
*       +ep^2*(12)
*       +ep^3*(32-14/3*z3)
*       +ep^4*(80-13/2*z4-56/3*z3)
*       +ep^5*(192-62/5*z5-26*z4-56*z3)
*       +ep^6*(448-61/3*z6-248/5*z5-78*z4-448/3*z3+98/9*z3^2)
*       +ep^7*(1024-254/7*z7-244/3*z6-744/5*z5-208*z4-1120/3*z3+91/3*z3*
*          z4+392/9*z3^2));
* id	GschemeConstants(0,0) = epp*(1
*       +ep*(2)
*       +ep^2*(4)
*       +ep^3*(8-7/3*z3)
*       +ep^4*(16-13/4*z4-14/3*z3)
*       +ep^5*(32-31/5*z5-13/2*z4-28/3*z3)
*       +ep^6*(64-61/6*z6-62/5*z5-13*z4-56/3*z3+49/18*z3^2)
*       +ep^7*(128-127/7*z7-61/3*z6-124/5*z5-26*z4-112/3*z3+91/12*z3*z4+
*          49/9*z3^2));


* id	epp = 1/ep;
#endprocedure

#procedure GammaArgToOne
*
* redcut
*

* 
* Reduce all Gam[a,b] to Gam[1,x] 
*         
        id Gam(x?,y?) * iGam(x?,y?) = 1;
        
        id Gam(x?pos_,0)  = fac_(x-1);
        id iGam(x?pos_,0) = 1/fac_(x-1);
        id Gam(0,y?)      = rat(1,ep*y) * Gam(1,y);
        id iGam(0,y?)     = rat(ep*y,1) * iGam(1,y);
        
        .sort
        
        repeat;
                id Gam(x?pos_,y?)         = rat(x-1+ep*y,1) * Gam(x-1,y);
                id Gam(1,y?) * iGam(1,y?) = 1;
        endrepeat;
        
        .sort
        
        repeat; 
                id iGam(x?pos_,y?)        = rat(1,x-1+ep*y) * iGam(x-1,y);
                id Gam(1,y?) * iGam(1,y?) = 1;
        endrepeat;
        
        .sort
        
        repeat;
                id Gam(x?neg0_,y?)        = rat(1,x+ep*y) * Gam(x+1,y);
                id Gam(1,y?) * iGam(1,y?) = 1;
        endrepeat;
        
        .sort

        
        repeat;
                id iGam(x?neg0_,y?)       = rat(x+ep*y,1) * iGam(x+1,y);
                id Gam(1,y?) * iGam(1,y?) = 1;
        endrepeat;
        
        .sort:All Gam and iGam reduced;
#endprocedure

* #procedure subMasters()

* * Two-loop fully massive tadpole        

* * MMM(1,1,1)        
* * id miT1 = M^2*(
* *         -21/2 - 3/(2*ep^2) - 9/(2*ep) + (27*S2)/2 - (3*z2)/2 
* *         + ep * T1ep + ep^2 * T1ep2
* *         );        

* * MMMMMM(1,1,1,1,1,1)        
* id miD6 = 2*z3/ep + D6 + ep*D6ep;

* * MMMMM0(1,1,1,1,1,1)        
* id miD5 = 2*z3/ep + D5 + ep*D5ep;

* * 00MMMM(1,1,1,1,1,1)        
* id miD4 = 2*z3/ep + D4 + ep*D4ep;

* * MM0000(1,1,1,1,1,1)        
* id miDN = 2*z3/ep + DN + ep*DNep;

* * MMM000(1,1,1,1,1,0)/M^2        
* *** the coefficients of 1/ep^3 and 1/ep^2 are determined
* *** from the cancellation of poles 
* id miE3 = -2/3/ep^3-11/3/ep^2 + (-14 + (27*S2)/2 - 2*z2)/ep + E3 + ep*E3ep;

* #endprocedure

#procedure cutep(x)
multiply, 1/ep^('x');
id ep=0;
multiply, ep^('x');
#endprocedure



#procedure subvalues
        
* Expansion from original MATAD        
        repeat  id,once    iGam(1,y?) = 
        1 - ep^2*y^2*z2/2 + ep^3*y^3*z3/3 + ep^4*y^4*(z2^2 - 2*z4)/8 +
        ep^5*y^5*(-5*z2*z3 + 6*z5)/30 +
        ep^6*y^6*(-3*z2^3 + 8*z3^2 + 18*z2*z4 - 24*z6)/144 +
        ep^7*iGamtrunc;
        
        
        repeat  id,once     Gam(1,y?) = 
        1 + ep^2*y^2*z2/2 - ep^3*y^3*z3/3 + ep^4*y^4*(z2^2 + 2*z4)/8 
        - ep^5*y^5*(5*z2*z3 + 6*z5)/30 +
        ep^6*y^6*(3*z2^3 + 8*z3^2 + 18*z2*z4 + 24*z6)/144 + 
        ep^7*Gamtrunc;
        


* Two-loop fully massive tadpole        
* MMM(1,1,1)/M^2        
        id miT1 = 
        -21/2 - 3/(2*ep^2) - 9/(2*ep) + (27*S2)/2 - (3*z2)/2 
        + ep * T1ep + ep^2 * T1ep2 + ep^3*miT1trunc;        
        
* MMMMMM(1,1,1,1,1,1)        
        id miD6 = 2*z3/ep + D6 + ep*D6ep + ep^3*miD6trunc;
        
* MMMMM0(1,1,1,1,1,1)        
        id miD5 = 2*z3/ep + D5 + ep*D5ep + ep^3*miD5trunc;
        
* 00MMMM(1,1,1,1,1,1)        
        id miD4 = 2*z3/ep + D4 + ep*D4ep + ep^3*miD4trunc;
        
* MM0000(1,1,1,1,1,1)        
        id miDN = 2*z3/ep + DN + ep*DNep + ep^3*miDNtrunc;
        
* MMM000(1,1,1,1,1,0)/M^2        
*** the coefficients of 1/ep^3 and 1/ep^2 are determined
*** from the cancellation of poles 
        id miE3 = -2/3/ep^3-11/3/ep^2 + (-14 + (27*S2)/2 - 2*z2)/ep + E3 + ep*E3ep + ep^3*miE3trunc;
        
#endprocedure


*--#[ expansion :
*
#procedure expansion(maxeppow)
*
*	Expands the PolyRatFun to sufficient powers in ep.
*
S DUMMYSYMBOL;
.sort:expansion-1;
PolyRatFun;
* 
* We introduce DUMMYSYMBOL for proper terms order in SplitArg
* First term has smallest power of ep        
*         
id	den(x?) = den(DUMMYSYMBOL*x);
id	rat(x1?,x2?) = num(x1)*den(DUMMYSYMBOL*x2);

SplitArg,den;
Multiply replace_(DUMMYSYMBOL,1);

* id	den(?a,x1?) = den(x1,?a);
repeat id den(x1?,x2?,x3?,?a) = den(x1,x2+x3,?a);
id	den(x1?,x2?) = den(1,x2/x1)/x1;
id	den(x1?) = 1/x1;

* Print+s;
* .end
.sort:expansion-2;
id	num(x1?) = x1;
if ( count(ep,1) > `maxeppow' ) discard;
repeat;
	id den(1,x?) = 1-x*den(1,x);
	if ( count(ep,1) > `maxeppow' ) discard;
endrepeat;
.sort:expansion-3;
Symbol ep(:`maxeppow');
#call subvalues
*
#endprocedure
*
*--#] expansion : 
