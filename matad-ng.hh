S n,ep;
dimension n;

* By default we use reduction with tables for BN topology

#ifndef `REDBNTAB'
        #define REDBNTAB
#endif

V q, Q;

S M,m;

V P,p,p1,...,p19;
S e1,...,e19;
V k,l,v,v1,v2,v3,v4,v5;
S s1m,s2m,s3m,s4m,s5m,s6m,s7m,s8m,s9m,s10m,s11m,s12m,s13m,s14m,s15m,s16m,s17m,s18m,s19m;


I i1,...,i9,j1,...,j9;
I MU,NU,al,be,ro,si,tau,mu,nu,la,ka;




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
S miT1,miD6,miD5,miD4,miDN,miDM,miE3;
S miBN,miBN1x11,miBN1x00;
* And its truncation flags
S miT1trunc,miD6trunc,miD5trunc,miD4trunc,miDNtrunc,miDMtrunc,miE3trunc;
S miBNtrunc,miBN1x11trunc,miBN1x00trunc;
S iGamtrunc,Gamtrunc;

set trunc:miT1trunc,miD6trunc,miD5trunc,miD4trunc,miDNtrunc,miE3trunc,  iGamtrunc,Gamtrunc;

* Mass distributions
* Two-loop 
Symbols   [000],[M00],[0M0],[00M],[MM0],[M0M],[0MM],[MMM];
* Three-loop 
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

.global



#procedure ToAuxTopo

*
*
*       AUX topo is T2
*
*       /----1->--\
*      |           |
*       ---<-3-----   
*      |           |
*       \--<-2----/        
*
*         
        id tad2l([000],n1?,n2?,n3?) = 1/p1.p1^n1/p2.p2^n2/p3.p3^n3;
        id tad2l([M00],n1?,n2?,n3?) = 1*s1m^n1/p2.p2^n2/p3.p3^n3;
        id tad2l([0M0],n1?,n2?,n3?) = 1/p1.p1^n1*s2m^n2/p3.p3^n3;
        id tad2l([00M],n1?,n2?,n3?) = 1/p1.p1^n1/p2.p2^n2*s3m^n3;
        id tad2l([MM0],n1?,n2?,n3?) = 1*s1m^n1*s2m^n2/p3.p3^n3;
        id tad2l([M0M],n1?,n2?,n3?) = 1*s1m^n1/p2.p2^n2*s3m^n3;
        id tad2l([0MM],n1?,n2?,n3?) = 1/p1.p1^n1*s2m^n2*s3m^n3;
        id tad2l([MMM],n1?,n2?,n3?) = 1*s1m^n1*s2m^n2*s3m^n3;                
        

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

        id po(1,?a) = 1;
        id poinv(1,?a) = 1;
        id po(x1?pos_,0) = fac_(x1-1);
        id poinv(x1?pos_,0) = 1/(fac_(x1-1));
        
        id po(n?,x?)                = Pochhammer(n,x*ep)*den(x*ep);
        id poinv(n?,x?)             = PochhammerINV(n,x*ep)*num(x*ep);
        
        id nom(x1?,x2?)             = num(x1+x2*ep);
        id deno(x1?,x2?)            = den(x1+x2*ep);

        id nom(x?,y?,z?)            = x + y*num(ep)*num(1+ep*z/y);

        id num(x?)                  = rat(x,1);
        id den(x?)                  = rat(1,x);        
        
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
        #call subSimple
        #call GammaArgToOne

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

*         #if `TYPE'==bn3
*                 Print+s;
*                 .end
                
*         #endif                
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
        
        id n = num(4-2*ep);
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

        if(count(intbn1,1));        
        id BN1(0,0,1,1,1,1) = +  M^4*miBN1x00*int0/intbn1;
        
        id BN1(1,1,1,1,1,1) = + miBN1x11*int0/intbn1;
        endif;
        .sort
        
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
        
        #do type = {d6|d5|d4|dm|dn|e4|e3}
                
                
                #message Recursion of type `type'
                
                #call top`type'
                Print+s;        
                .sort        
                
        #enddo
        
        #do type = {bn|bn1|bn2|bn3|bm|bm1|bm2}
        
                #message Recursion of type `type'
                multiply replace_(s1m,x1,s2m,x2,s3m,x3,s4m,x4,s5m,x5,s6m,x6,s7m,x7,s8m,x8);
                
                #call top`type'
                #call bnm2m(`type')
                Print+s;        
                .sort        
                
        #enddo
        
        #do type = {m1|m2|m3|m4|m5|t1|n1}
                #message Recursion of type `type'
                multiply replace_(s1m,x1,s2m,x2,s3m,x3,s4m,x4,s5m,x5,s6m,x6,s7m,x7,s8m,x8);
                
                #call top`type'
        #enddo

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
#procedure topd6

* treat the scalar products


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
        id s1m*s2m*s3m/p4.p4/p5.p5/p6.p6 = miDM*int0/intdm;
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
id s1m*s2m/p3.p3/p4.p4/p5.p5/p6.p6 = miDN*int0/intdn;
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


****************************************************************************
* 
* New tables for BN(1,1,1,1) BN(1,1,1,2) BN(1,1,2,2) BN(1,2,2,2) BN(2,2,2,2)
* and it's dalambertian upto dala^11
* 
#procedure BNd0Exact
id,only intbn*dala^0*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^4*miBN*rat(1,1) );


id,only intbn*dala^0*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^2*miBN*rat(3*ep - 2,4
         ) );


id,only intbn*dala^0*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( miBN*rat(18*ep^3 - 27*
         ep^2 + 13*ep - 2,32*ep) + Gam(1,1)^3*rat(1,8*ep^4) );


id,only intbn*dala^0*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-2*miBN*rat(54*ep^4
          - 135*ep^3 + 120*ep^2 - 45*ep + 6,128*ep) + Gam(1,1)^3*M^-2*rat(11*
         ep - 3,32*ep^4) );


id,only intbn*dala^0*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-4*miBN*rat(324*ep^6
          - 864*ep^5 + 315*ep^4 + 312*ep^3 - 147*ep^2 - 24*ep + 12,1024*ep^2
          + 1024*ep) + Gam(1,1)^3*M^-4*rat(162*ep^3 + 67*ep^2 - 11*ep - 6,256*
         ep^5 + 256*ep^4) );
#endprocedure

#procedure BNdExact
id,only intbn*dala^1*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^2*miBN*rat(6*ep^3 - 
         25*ep^2 + 23*ep - 6,8*ep) + Gam(1,1)^3*M^2*rat(3,2*ep^4) );


id,only intbn*dala^1*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( miBN*rat(18*ep^4 - 81*
         ep^3 + 94*ep^2 - 41*ep + 6,32*ep) + Gam(1,1)^3*rat(9*ep - 3,8*ep^4) )
         ;


id,only intbn*dala^1*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-2*miBN*rat(108*ep^6
          - 432*ep^5 + 321*ep^4 + 144*ep^3 - 249*ep^2 + 96*ep - 12,256*ep^2 + 
         256*ep) + Gam(1,1)^3*M^-2*rat(54*ep^3 - 7*ep^2 - 25*ep + 6,64*ep^5 + 
         64*ep^4) );


id,only intbn*dala^1*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-4*miBN*rat(324*ep^7
          - 1512*ep^6 + 2043*ep^5 - 318*ep^4 - 771*ep^3 + 270*ep^2 + 60*ep - 
         24,1024*ep^2 + 1024*ep) + Gam(1,1)^3*M^-4*rat(162*ep^4 - 257*ep^3 - 
         145*ep^2 + 16*ep + 12,256*ep^5 + 256*ep^4) );


id,only intbn*dala^1*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-6*miBN*rat(1944*ep^9
          - 6804*ep^8 - 594*ep^7 + 32697*ep^6 - 417*ep^5 - 23205*ep^4 + 159*
         ep^3 + 5208*ep^2 - 12*ep - 336,8192*ep^3 + 24576*ep^2 + 16384*ep) + 
         Gam(1,1)^3*M^-6*rat(972*ep^6 - 1944*ep^5 - 9947*ep^4 - 8084*ep^3 - 
         1111*ep^2 + 650*ep + 168,2048*ep^6 + 6144*ep^5 + 4096*ep^4) );


id,only intbn*dala^2*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( miBN*rat(36*ep^6 - 288*
         ep^5 + 755*ep^4 - 416*ep^3 - 187*ep^2 + 192*ep - 36,64*ep^2 + 64*ep)
          + Gam(1,1)^3*rat(18*ep^3 - 117*ep^2 - 27*ep + 18,16*ep^5 + 16*ep^4)
          );


id,only intbn*dala^2*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-2*miBN*rat(108*ep^6
          - 864*ep^5 + 2265*ep^4 - 1248*ep^3 - 561*ep^2 + 576*ep - 108,256*ep
          + 256) + Gam(1,1)^3*M^-2*rat(54*ep^3 - 351*ep^2 - 81*ep + 54,64*ep^4
          + 64*ep^3) );


id,only intbn*dala^2*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-4*miBN*rat(648*ep^9
          - 3996*ep^8 + 4410*ep^7 + 16131*ep^6 - 20235*ep^5 - 231*ep^4 + 6837*
         ep^3 - 1464*ep^2 - 516*ep + 144,2048*ep^3 + 6144*ep^2 + 4096*ep) + 
         Gam(1,1)^3*M^-4*rat(324*ep^6 - 1512*ep^5 - 3673*ep^4 + 820*ep^3 + 
         1091*ep^2 - 18*ep - 72,512*ep^6 + 1536*ep^5 + 1024*ep^4) );


id,only intbn*dala^2*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-6*miBN*rat(1944*
         ep^10 - 12636*ep^9 + 19818*ep^8 + 34479*ep^7 - 98508*ep^6 - 21954*
         ep^5 + 69774*ep^4 + 4731*ep^3 - 15636*ep^2 - 300*ep + 1008,8192*ep^3
          + 24576*ep^2 + 16384*ep) + Gam(1,1)^3*M^-6*rat(972*ep^7 - 4860*ep^6
          - 4115*ep^5 + 21757*ep^4 + 23141*ep^3 + 3983*ep^2 - 1782*ep - 504,
         2048*ep^6 + 6144*ep^5 + 4096*ep^4) );


id,only intbn*dala^2*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-8*miBN*rat(11664*
         ep^12 - 46656*ep^11 - 70632*ep^10 + 644112*ep^9 - 493767*ep^8 - 
         4801392*ep^7 - 2931642*ep^6 + 3581892*ep^5 + 2604141*ep^4 - 817812*
         ep^3 - 625716*ep^2 + 53136*ep + 41472,65536*ep^4 + 393216*ep^3 + 
         720896*ep^2 + 393216*ep) + Gam(1,1)^3*M^-8*rat(5832*ep^9 - 14580*ep^8
          - 48438*ep^7 + 556593*ep^6 + 2043681*ep^5 + 2488773*ep^4 + 1138125*
         ep^3 + 29574*ep^2 - 106056*ep - 20736,16384*ep^7 + 98304*ep^6 + 
         180224*ep^5 + 98304*ep^4) );


id,only intbn*dala^3*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-2*miBN*rat(216*ep^9
          - 2484*ep^8 + 9822*ep^7 - 8415*ep^6 - 38049*ep^5 + 32499*ep^4 + 8175
         *ep^3 - 9288*ep^2 + 180*ep + 432,512*ep^3 + 1536*ep^2 + 1024*ep) + 
         Gam(1,1)^3*M^-2*rat(108*ep^6 - 1080*ep^5 + 3453*ep^4 + 10452*ep^3 + 
         2409*ep^2 - 918*ep - 216,128*ep^6 + 384*ep^5 + 256*ep^4) );


id,only intbn*dala^3*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-4*miBN*rat(648*ep^10
          - 7236*ep^9 + 26982*ep^8 - 15423*ep^7 - 122562*ep^6 + 59448*ep^5 + 
         57024*ep^4 - 19689*ep^3 - 8748*ep^2 + 1476*ep + 432,2048*ep^3 + 6144*
         ep^2 + 4096*ep) + Gam(1,1)^3*M^-4*rat(324*ep^7 - 3132*ep^6 + 9279*
         ep^5 + 34809*ep^4 + 17679*ep^3 - 345*ep^2 - 1566*ep - 216,512*ep^6 + 
         1536*ep^5 + 1024*ep^4) );


id,only intbn*dala^3*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-6*miBN*rat(3888*
         ep^12 - 31104*ep^11 + 30888*ep^10 + 359424*ep^9 - 916053*ep^8 - 
         1639896*ep^7 + 1562610*ep^6 + 1162668*ep^5 - 844377*ep^4 - 261060*
         ep^3 + 169092*ep^2 + 16848*ep - 10368,16384*ep^4 + 98304*ep^3 + 
         180224*ep^2 + 98304*ep) + Gam(1,1)^3*M^-6*rat(1944*ep^9 - 12636*ep^8
          - 594*ep^7 + 335091*ep^6 + 605931*ep^5 + 122295*ep^4 - 207873*ep^3
          - 71478*ep^2 + 11448*ep + 5184,4096*ep^7 + 24576*ep^6 + 45056*ep^5
          + 24576*ep^4) );


id,only intbn*dala^3*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-8*miBN*rat(11664*
         ep^13 - 93312*ep^12 + 115992*ep^11 + 926640*ep^10 - 3070215*ep^9 - 
         2826324*ep^8 + 16273926*ep^7 + 15308460*ep^6 - 11723427*ep^5 - 
         11234376*ep^4 + 2645532*ep^3 + 2556000*ep^2 - 171072*ep - 165888,
         65536*ep^4 + 393216*ep^3 + 720896*ep^2 + 393216*ep) + Gam(1,1)^3*M^-8
         *rat(5832*ep^10 - 37908*ep^9 + 9882*ep^8 + 750345*ep^7 - 182691*ep^6
          - 5685951*ep^5 - 8816967*ep^4 - 4522926*ep^3 - 224352*ep^2 + 403488*
         ep + 82944,16384*ep^7 + 98304*ep^6 + 180224*ep^5 + 98304*ep^4) );


id,only intbn*dala^3*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-10*miBN*rat(69984*
         ep^15 - 291600*ep^14 - 1321920*ep^13 + 8530920*ep^12 - 4930578*ep^11
          - 117134649*ep^10 + 256737339*ep^9 + 1783904976*ep^8 + 2257021044*
         ep^7 - 244562517*ep^6 - 1935136953*ep^5 - 598519602*ep^4 + 460273500*
         ep^3 + 193809816*ep^2 - 30386016*ep - 14079744,524288*ep^5 + 5242880*
         ep^4 + 18350080*ep^3 + 26214400*ep^2 + 12582912*ep) + Gam(1,1)^3*
         M^-10*rat(34992*ep^12 - 93312*ep^11 - 748440*ep^10 + 2617920*ep^9 - 
         21490989*ep^8 - 280174056*ep^7 - 949448058*ep^6 - 1500340416*ep^5 - 
         1182911037*ep^4 - 388001544*ep^3 + 22933836*ep^2 + 42179184*ep + 
         7039872,131072*ep^8 + 1310720*ep^7 + 4587520*ep^6 + 6553600*ep^5 + 
         3145728*ep^4) );


id,only intbn*dala^4*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-4*miBN*rat(1296*
         ep^12 - 19008*ep^11 + 92952*ep^10 - 54480*ep^9 - 1095231*ep^8 + 
         2713416*ep^7 + 5853750*ep^6 - 2700012*ep^5 - 3870219*ep^4 + 698148*
         ep^3 + 825804*ep^2 - 47952*ep - 51840,4096*ep^4 + 24576*ep^3 + 45056*
         ep^2 + 24576*ep) + Gam(1,1)^3*M^-4*rat(648*ep^9 - 8532*ep^8 + 34650*
         ep^7 + 4809*ep^6 - 1138239*ep^5 - 2429883*ep^4 - 1425075*ep^3 - 94194
         *ep^2 + 123336*ep + 25920,1024*ep^7 + 6144*ep^6 + 11264*ep^5 + 6144*
         ep^4) );


id,only intbn*dala^4*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-6*miBN*rat(3888*
         ep^13 - 54432*ep^12 + 240840*ep^11 + 22464*ep^10 - 3394653*ep^9 + 
         5949786*ep^8 + 22988082*ep^7 + 3607464*ep^6 - 17010681*ep^5 - 5645994
         *ep^4 + 3873708*ep^3 + 1507752*ep^2 - 251424*ep - 103680,16384*ep^4
          + 98304*ep^3 + 180224*ep^2 + 98304*ep) + Gam(1,1)^3*M^-6*rat(1944*
         ep^10 - 24300*ep^9 + 86886*ep^8 + 83727*ep^7 - 3405099*ep^6 - 9566127
         *ep^5 - 9134991*ep^4 - 3132732*ep^3 + 181620*ep^2 + 324432*ep + 51840
         ,4096*ep^7 + 24576*ep^6 + 45056*ep^5 + 24576*ep^4) );


id,only intbn*dala^4*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-8*miBN*rat(23328*
         ep^15 - 221616*ep^14 + 57024*ep^13 + 5587704*ep^12 - 16476678*ep^11
          - 45832563*ep^10 + 246722409*ep^9 + 541636272*ep^8 - 105304068*ep^7
          - 609631911*ep^6 - 25012755*ep^5 + 251607834*ep^4 + 13284468*ep^3 - 
         42627960*ep^2 - 1060128*ep + 2384640,131072*ep^5 + 1310720*ep^4 + 
         4587520*ep^3 + 6553600*ep^2 + 3145728*ep) + Gam(1,1)^3*M^-8*rat(11664
         *ep^12 - 93312*ep^11 - 93960*ep^10 + 2384640*ep^9 - 15588975*ep^8 - 
         123457368*ep^7 - 248037726*ep^6 - 163465248*ep^5 + 30131457*ep^4 + 
         66037416*ep^3 + 12913092*ep^2 - 4040496*ep - 1192320,32768*ep^8 + 
         327680*ep^7 + 1146880*ep^6 + 1638400*ep^5 + 786432*ep^4) );


id,only intbn*dala^4*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-10*miBN*rat(69984*
         ep^16 - 641520*ep^15 + 136080*ep^14 + 15140520*ep^13 - 47585178*ep^12
          - 92481759*ep^11 + 842410584*ep^10 + 500218281*ep^9 - 6662503836*
         ep^8 - 11529667737*ep^7 - 712324368*ep^6 + 9077165163*ep^5 + 
         3452871510*ep^4 - 2107557684*ep^3 - 999435096*ep^2 + 137850336*ep + 
         70398720,524288*ep^5 + 5242880*ep^4 + 18350080*ep^3 + 26214400*ep^2
          + 12582912*ep) + Gam(1,1)^3*M^-10*rat(34992*ep^13 - 268272*ep^12 - 
         281880*ep^11 + 6360120*ep^10 - 34580589*ep^9 - 172719111*ep^8 + 
         451422222*ep^7 + 3246899874*ep^6 + 6318791043*ep^5 + 5526553641*ep^4
          + 1962941556*ep^3 - 72489996*ep^2 - 203856048*ep - 35199360,131072*
         ep^8 + 1310720*ep^7 + 4587520*ep^6 + 6553600*ep^5 + 3145728*ep^4) );


id,only intbn*dala^4*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-12*miBN*rat(419904*
         ep^17 - 2519424*ep^16 - 12212208*ep^15 + 113724000*ep^14 - 151059492*
         ep^13 - 1683245448*ep^12 + 10938699099*ep^11 + 17484608970*ep^10 - 
         227981745153*ep^9 - 795700856682*ep^8 - 785497422699*ep^7 + 
         219909241710*ep^6 + 735263239041*ep^5 + 185516095554*ep^4 - 
         178814773692*ep^3 - 68171800680*ep^2 + 11903047200*ep + 5103648000,
         4194304*ep^5 + 54525952*ep^4 + 247463936*ep^3 + 448790528*ep^2 + 
         251658240*ep) + Gam(1,1)^3*M^-12*rat(209952*ep^14 - 944784*ep^13 - 
         7208352*ep^12 + 42322824*ep^11 + 1391418*ep^10 + 1138887039*ep^9 + 
         30107046219*ep^8 + 171318360246*ep^7 + 448167167376*ep^6 + 
         621885198831*ep^5 + 453001390719*ep^4 + 138754523316*ep^3 - 
         11056733460*ep^2 - 15733515600*ep - 2551824000,1048576*ep^8 + 
         13631488*ep^7 + 61865984*ep^6 + 112197632*ep^5 + 62914560*ep^4) );


id,only intbn*dala^5*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-6*miBN*rat(7776*
         ep^14 - 143856*ep^13 + 878256*ep^12 - 443016*ep^11 - 19770042*ep^10
          + 83648697*ep^9 + 91014738*ep^8 - 1117956018*ep^7 - 1450097994*ep^6
          + 862161549*ep^5 + 1130693706*ep^4 - 198865116*ep^3 - 261816840*ep^2
          + 12972960*ep + 17107200,32768*ep^4 + 294912*ep^3 + 851968*ep^2 + 
         786432*ep) + Gam(1,1)^3*M^-6*rat(3888*ep^11 - 66096*ep^10 + 345816*
         ep^9 + 155304*ep^8 - 8149677*ep^7 + 82446525*ep^6 + 488308833*ep^5 + 
         793646535*ep^4 + 446481468*ep^3 + 31199580*ep^2 - 39275280*ep - 
         8553600,8192*ep^7 + 73728*ep^6 + 212992*ep^5 + 196608*ep^4) );


id,only intbn*dala^5*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-8*miBN*rat(23328*
         ep^15 - 408240*ep^14 + 2203200*ep^13 + 1305720*ep^12 - 60639174*ep^11
          + 191635965*ep^10 + 523990305*ep^9 - 3080823840*ep^8 - 7704162036*
         ep^7 - 1763809335*ep^6 + 5978565765*ep^5 + 2795485770*ep^4 - 
         1382045868*ep^3 - 746531640*ep^2 + 90240480*ep + 51321600,131072*ep^4
          + 1179648*ep^3 + 3407872*ep^2 + 3145728*ep) + Gam(1,1)^3*M^-8*rat(
         11664*ep^12 - 186624*ep^11 + 839160*ep^10 + 1503360*ep^9 - 23983119*
         ep^8 + 222890544*ep^7 + 1712266074*ep^6 + 3845866104*ep^5 + 
         3720384009*ep^4 + 1433043144*ep^3 - 24227100*ep^2 - 143486640*ep - 
         25660800,32768*ep^7 + 294912*ep^6 + 851968*ep^5 + 786432*ep^4) );


id,only intbn*dala^5*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-10*miBN*rat(139968*
         ep^18 - 1492992*ep^17 - 1862352*ep^16 + 70720128*ep^15 - 185482764*
         ep^14 - 1036687104*ep^13 + 6929504721*ep^12 + 8053978536*ep^11 - 
         127167791751*ep^10 - 361735949916*ep^9 - 194403886461*ep^8 + 
         327780306144*ep^7 + 318153329847*ep^6 - 99168439788*ep^5 - 
         138018192408*ep^4 + 11895864192*ep^3 + 23310817200*ep^2 - 492091200*
         ep - 1290816000,1048576*ep^6 + 15728640*ep^5 + 89128960*ep^4 + 
         235929600*ep^3 + 287309824*ep^2 + 125829120*ep) + Gam(1,1)^3*M^-10*
         rat(69984*ep^15 - 641520*ep^14 - 1788480*ep^13 + 30945240*ep^12 - 
         37330578*ep^11 + 417623145*ep^10 + 15403234179*ep^9 + 82673128668*
         ep^8 + 190214482836*ep^7 + 205362639021*ep^6 + 70650478047*ep^5 - 
         46357028562*ep^4 - 41078799396*ep^3 - 5064913800*ep^2 + 2720109600*ep
          + 645408000,262144*ep^9 + 3932160*ep^8 + 22282240*ep^7 + 58982400*
         ep^6 + 71827456*ep^5 + 31457280*ep^4) );


id,only intbn*dala^5*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-12*miBN*rat(419904*
         ep^18 - 5038848*ep^17 + 2904336*ep^16 + 186997248*ep^15 - 833403492*
         ep^14 - 776888496*ep^13 + 21038171787*ep^12 - 48147585624*ep^11 - 
         332889398973*ep^10 + 572189614236*ep^9 + 3988707717393*ep^8 + 
         4932893777904*ep^7 - 584192211219*ep^6 - 4226063338692*ep^5 - 
         1291911347016*ep^4 + 1004716841472*ep^3 + 420933851280*ep^2 - 
         66314635200*ep - 30621888000,4194304*ep^5 + 54525952*ep^4 + 247463936
         *ep^3 + 448790528*ep^2 + 251658240*ep) + Gam(1,1)^3*M^-12*rat(209952*
         ep^15 - 2204496*ep^14 - 1539648*ep^13 + 85572936*ep^12 - 252545526*
         ep^11 + 1130538531*ep^10 + 23273723985*ep^9 - 9323917068*ep^8 - 
         579742994100*ep^7 - 2067117805425*ep^6 - 3278309802267*ep^5 - 
         2579253820998*ep^4 - 843583873356*ep^3 + 50606885160*ep^2 + 
         91849269600*ep + 15310944000,1048576*ep^8 + 13631488*ep^7 + 61865984*
         ep^6 + 112197632*ep^5 + 62914560*ep^4) );


id,only intbn*dala^5*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-14*miBN*rat(2519424*
         ep^21 - 8817984*ep^20 - 187907040*ep^19 + 771258672*ep^18 + 
         3351784536*ep^17 - 26789487372*ep^16 + 103817077542*ep^15 + 
         554312161809*ep^14 - 7502377161357*ep^13 - 21096357822135*ep^12 + 
         212264679038859*ep^11 + 1407261578904711*ep^10 + 3530608360711785*
         ep^9 + 4020728271593247*ep^8 + 769249291495491*ep^7 - 
         2785063035366612*ep^6 - 2294416442972088*ep^5 + 67526202621744*ep^4
          + 652718793105648*ep^3 + 143071014121920*ep^2 - 45795094156800*ep - 
         13546662912000,33554432*ep^7 + 704643072*ep^6 + 5872025600*ep^5 + 
         24662507520*ep^4 + 54492397568*ep^3 + 59190018048*ep^2 + 24159191040*
         ep) + Gam(1,1)^3*M^-14*rat(1259712*ep^18 - 2519424*ep^17 - 95843088*
         ep^16 + 224228736*ep^15 + 1979325396*ep^14 - 9893002392*ep^13 - 
         153420682707*ep^12 - 3988501251732*ep^11 - 47776751165313*ep^10 - 
         289752074018838*ep^9 - 1033706547769185*ep^8 - 2323452765971376*ep^7
          - 3353799816489411*ep^6 - 3038483193155622*ep^5 - 1576463447412396*
         ep^4 - 322614213750504*ep^3 + 75505073646240*ep^2 + 48861984326400*ep
          + 6773331456000,8388608*ep^10 + 176160768*ep^9 + 1468006400*ep^8 + 
         6165626880*ep^7 + 13623099392*ep^6 + 14797504512*ep^5 + 6039797760*
         ep^4) );


id,only intbn*dala^6*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-8*miBN*rat(46656*
         ep^18 - 933120*ep^17 + 5071248*ep^16 + 16568064*ep^15 - 287305668*
         ep^14 + 877909680*ep^13 + 4160205891*ep^12 - 36038356632*ep^11 - 
         37078448517*ep^10 + 683291475900*ep^9 + 2135362771449*ep^8 + 
         1885778793072*ep^7 - 725209849611*ep^6 - 1825601593860*ep^5 - 
         391303622088*ep^4 + 447499610496*ep^3 + 158262938640*ep^2 - 
         29874225600*ep - 12083904000,262144*ep^6 + 3932160*ep^5 + 22282240*
         ep^4 + 58982400*ep^3 + 71827456*ep^2 + 31457280*ep) + Gam(1,1)^3*M^-8
         *rat(23328*ep^15 - 431568*ep^14 + 1923264*ep^13 + 10264968*ep^12 - 
         119083590*ep^11 + 213280443*ep^10 - 3433673511*ep^9 - 89473458060*
         ep^8 - 482552416116*ep^7 - 1203088078521*ep^6 - 1598373481779*ep^5 - 
         1118723407254*ep^4 - 327156667308*ep^3 + 30994543080*ep^2 + 
         38097928800*ep + 6041952000,65536*ep^9 + 983040*ep^8 + 5570560*ep^7
          + 14745600*ep^6 + 17956864*ep^5 + 7864320*ep^4) );


id,only intbn*dala^6*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-10*miBN*rat(139968*
         ep^19 - 2612736*ep^18 + 11481264*ep^17 + 69989184*ep^16 - 795644748*
         ep^15 + 1484506368*ep^14 + 15992256393*ep^13 - 91474246332*ep^12 - 
         255388772079*ep^11 + 1901560633632*ep^10 + 9139254217947*ep^9 + 
         14198787465012*ep^8 + 5367485623455*ep^7 - 8377644180024*ep^6 - 
         8476317241704*ep^5 - 222715656864*ep^4 + 2264787257904*ep^3 + 
         543429077760*ep^2 - 155748614400*ep - 48335616000,1048576*ep^6 + 
         15728640*ep^5 + 89128960*ep^4 + 235929600*ep^3 + 287309824*ep^2 + 
         125829120*ep) + Gam(1,1)^3*M^-10*rat(69984*ep^16 - 1201392*ep^15 + 
         4043520*ep^14 + 38487960*ep^13 - 316190898*ep^12 + 163506969*ep^11 - 
         9447898761*ep^10 - 282155068224*ep^9 - 1805551080588*ep^8 - 
         5539473900027*ep^7 - 9607472759421*ep^6 - 9749664148878*ep^5 - 
         5456363630940*ep^4 - 1215643039992*ep^3 + 238271958720*ep^2 + 
         170517571200*ep + 24167808000,262144*ep^9 + 3932160*ep^8 + 22282240*
         ep^7 + 58982400*ep^6 + 71827456*ep^5 + 31457280*ep^4) );


id,only intbn*dala^6*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-12*miBN*rat(839808*
         ep^21 - 9657792*ep^20 - 35761824*ep^19 + 774967824*ep^18 - 1278990648
         *ep^17 - 20088653508*ep^16 + 117808083810*ep^15 + 107615169963*ep^14
          - 4274154414639*ep^13 + 417673086291*ep^12 + 113454038693529*ep^11
          + 395942363470557*ep^10 + 477528159364515*ep^9 - 47652452733915*ep^8
          - 548305558403967*ep^7 - 254582751663324*ep^6 + 195887862992664*ep^5
          + 135217015505424*ep^4 - 28151201328048*ep^3 - 24150737891520*ep^2
          + 1386980236800*ep + 1363848192000,8388608*ep^7 + 176160768*ep^6 + 
         1468006400*ep^5 + 6165626880*ep^4 + 13623099392*ep^3 + 14797504512*
         ep^2 + 6039797760*ep) + Gam(1,1)^3*M^-12*rat(419904*ep^18 - 4199040*
         ep^17 - 23549616*ep^16 + 341241984*ep^15 - 89053668*ep^14 - 
         9914662440*ep^13 - 53333020449*ep^12 - 1798100191980*ep^11 - 
         18972501116715*ep^10 - 90583230387762*ep^9 - 235497669047739*ep^8 - 
         346129564525920*ep^7 - 260154937695393*ep^6 - 39750417955170*ep^5 + 
         76897943454492*ep^4 + 46482762820104*ep^3 + 3450154318560*ep^2 - 
         3307532486400*ep - 681924096000,2097152*ep^10 + 44040192*ep^9 + 
         367001600*ep^8 + 1541406720*ep^7 + 3405774848*ep^6 + 3699376128*ep^5
          + 1509949440*ep^4) );


id,only intbn*dala^6*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-14*miBN*rat(2519424*
         ep^22 - 26453952*ep^21 - 126181152*ep^20 + 2086607952*ep^19 - 
         2047026168*ep^18 - 50251979124*ep^17 + 291343489146*ep^16 - 
         172407380985*ep^15 - 11382562294020*ep^14 + 31420282307364*ep^13 + 
         359939183793804*ep^12 - 78591174367302*ep^11 - 6320222691621192*ep^10
          - 20693530253389248*ep^9 - 27375848609657238*ep^8 - 8169808075835049
         *ep^7 + 17201024804594196*ep^6 + 16128441303426360*ep^5 + 
         180035374753440*ep^4 - 4425960537617616*ep^3 - 1047292193010240*ep^2
          + 307018996185600*ep + 94826640384000,33554432*ep^7 + 704643072*ep^6
          + 5872025600*ep^5 + 24662507520*ep^4 + 54492397568*ep^3 + 
         59190018048*ep^2 + 24159191040*ep) + Gam(1,1)^3*M^-14*rat(1259712*
         ep^19 - 11337408*ep^18 - 78207120*ep^17 + 895130352*ep^16 + 409724244
         *ep^15 - 23748280164*ep^14 - 84169665963*ep^13 - 2914556472783*ep^12
          - 19857242403189*ep^11 + 44685184138353*ep^10 + 994557970362681*ep^9
          + 4912493068412919*ep^8 + 12910369545310221*ep^7 + 20438115522270255
         *ep^6 + 19692918904676958*ep^5 + 10712629918136268*ep^4 + 
         2333804569899768*ep^3 - 479673531197280*ep^2 - 335260558828800*ep - 
         47413320192000,8388608*ep^10 + 176160768*ep^9 + 1468006400*ep^8 + 
         6165626880*ep^7 + 13623099392*ep^6 + 14797504512*ep^5 + 6039797760*
         ep^4) );


id,only intbn*dala^6*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-16*miBN*rat(15116544
         *ep^24 - 40310784*ep^23 - 1829101824*ep^22 + 5394926592*ep^21 + 
         67745141856*ep^20 - 281974307328*ep^19 + 539718955056*ep^18 + 
         6033692279232*ep^17 - 152342545570119*ep^16 - 7181691855048*ep^15 + 
         9097172697965436*ep^14 + 14943011674960152*ep^13 - 328950853687087794
         *ep^12 - 2301616403966272248*ep^11 - 7140307187117150604*ep^10 - 
         11919900242908426968*ep^9 - 9277396851562325199*ep^8 + 
         1408610648195587872*ep^7 + 8430966650367874296*ep^6 + 
         4887946510400738592*ep^5 - 827458459095546288*ep^4 - 
         1571605382270542464*ep^3 - 265760210153007360*ep^2 + 
         114201282034022400*ep + 29714897350656000,268435456*ep^8 + 7516192768
         *ep^7 + 86436216832*ep^6 + 526133493760*ep^5 + 1817039601664*ep^4 + 
         3525094408192*ep^3 + 3507914539008*ep^2 + 1352914698240*ep) + Gam(1,1
         )^3*M^-16*rat(7558272*ep^21 - 8817984*ep^20 - 916440480*ep^19 + 
         1226434608*ep^18 + 34933405896*ep^17 - 81398025900*ep^16 + 
         174370679154*ep^15 + 24413519344689*ep^14 + 604921089938085*ep^13 + 
         9694105263680391*ep^12 + 90039439355115537*ep^11 + 517053506711573421
         *ep^10 + 1951999971378659307*ep^9 + 5013106757265984417*ep^8 + 
         8852598075971018457*ep^7 + 10623551524854856542*ep^6 + 
         8302392182260351980*ep^5 + 3776555787160435416*ep^4 + 
         625304218958220672*ep^3 - 216008361397825920*ep^2 - 
         114054194272435200*ep - 14857448675328000,67108864*ep^11 + 1879048192
         *ep^10 + 21609054208*ep^9 + 131533373440*ep^8 + 454259900416*ep^7 + 
         881273602048*ep^6 + 876978634752*ep^5 + 338228674560*ep^4) );


id,only intbn*dala^7*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-10*miBN*rat(279936*
         ep^21 - 6205248*ep^20 + 30754080*ep^19 + 270681264*ep^18 - 3380688360
         *ep^17 + 7392219156*ep^16 + 86610022614*ep^15 - 744794946567*ep^14 - 
         133886923821*ep^13 + 27514741243185*ep^12 - 3184892281557*ep^11 - 
         738666634740993*ep^10 - 2777751399430647*ep^9 - 4146602087714121*ep^8
          - 1658132540490573*ep^7 + 2403228943597356*ep^6 + 2616032460193992*
         ep^5 + 142396871844528*ep^4 - 699416922626064*ep^3 - 180920010922560*
         ep^2 + 48114251942400*ep + 15692092416000,2097152*ep^7 + 44040192*
         ep^6 + 367001600*ep^5 + 1541406720*ep^4 + 3405774848*ep^3 + 
         3699376128*ep^2 + 1509949440*ep) + Gam(1,1)^3*M^-10*rat(139968*ep^18
          - 2892672*ep^17 + 11247984*ep^16 + 146333952*ep^15 - 1412913996*
         ep^14 + 1404752904*ep^13 + 44816186565*ep^12 + 365853074580*ep^11 + 
         12611898920367*ep^10 + 130154054570874*ep^9 + 633927266752311*ep^8 + 
         1753695964092024*ep^7 + 2938584791517309*ep^6 + 2982183884312922*ep^5
          + 1701736203888468*ep^4 + 396069350533272*ep^3 - 70411881748320*ep^2
          - 54133636435200*ep - 7846046208000,524288*ep^10 + 11010048*ep^9 + 
         91750400*ep^8 + 385351680*ep^7 + 851443712*ep^6 + 924844032*ep^5 + 
         377487360*ep^4) );


id,only intbn*dala^7*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-12*miBN*rat(839808*
         ep^22 - 17216064*ep^21 + 61236000*ep^20 + 965814192*ep^19 - 
         8788658760*ep^18 + 5273215668*ep^17 + 296791163622*ep^16 - 
         1801334726631*ep^15 - 4125635504298*ep^14 + 81874789110450*ep^13 + 
         128019029371254*ep^12 - 2231924365630764*ep^11 - 12026587371996906*
         ep^10 - 26328563260295598*ep^9 - 25707408060042324*ep^8 - 
         1080975871660797*ep^7 + 19864242098568756*ep^6 + 13507352916503544*
         ep^5 - 1386266408655552*ep^4 - 4039844645898000*ep^3 - 
         760257298785600*ep^2 + 287647536960000*ep + 78460462080000,8388608*
         ep^7 + 176160768*ep^6 + 1468006400*ep^5 + 6165626880*ep^4 + 
         13623099392*ep^3 + 14797504512*ep^2 + 6039797760*ep) + Gam(1,1)^3*
         M^-12*rat(419904*ep^19 - 7978176*ep^18 + 19280592*ep^17 + 495241776*
         ep^16 - 3507072228*ep^15 - 2850311268*ep^14 + 141472324215*ep^13 + 
         1321640156565*ep^12 + 39664962134001*ep^11 + 453521658314457*ep^10 + 
         2552552073111303*ep^9 + 8430724226037627*ep^8 + 17584234195012047*
         ep^7 + 23639475610525311*ep^6 + 20016128033230014*ep^5 + 
         9696889071042156*ep^4 + 1769111107421400*ep^3 - 514460318047200*ep^2
          - 294206320800000*ep - 39230231040000,2097152*ep^10 + 44040192*ep^9
          + 367001600*ep^8 + 1541406720*ep^7 + 3405774848*ep^6 + 3699376128*
         ep^5 + 1509949440*ep^4) );


id,only intbn*dala^7*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-14*miBN*rat(5038848*
         ep^24 - 60466176*ep^23 - 445098240*ep^22 + 7610340096*ep^21 + 
         119369376*ep^20 - 328430386368*ep^19 + 1385421347952*ep^18 + 
         3139077567888*ep^17 - 85631459993397*ep^16 + 232608496652064*ep^15 + 
         3933200436002964*ep^14 - 5843444427326520*ep^13 - 155215302656563350*
         ep^12 - 638958730961167704*ep^11 - 1167980732603948484*ep^10 - 
         758329620479678808*ep^9 + 637148351890024947*ep^8 + 
         1215972674955217368*ep^7 + 261657404658357768*ep^6 - 
         481494820509529056*ep^5 - 216726178630769424*ep^4 + 74006503761378816
         *ep^3 + 42528421906871040*ep^2 - 3820152110745600*ep - 
         2485715226624000,67108864*ep^8 + 1879048192*ep^7 + 21609054208*ep^6
          + 131533373440*ep^5 + 454259900416*ep^4 + 881273602048*ep^3 + 
         876978634752*ep^2 + 338228674560*ep) + Gam(1,1)^3*M^-14*rat(2519424*
         ep^21 - 26453952*ep^20 - 258450912*ep^19 + 3350099088*ep^18 + 
         5154432408*ep^17 - 151744560516*ep^16 + 447887488230*ep^15 + 
         12812115487635*ep^14 + 267043587157143*ep^13 + 3982872325196637*ep^12
          + 31773702688472619*ep^11 + 147947222729676303*ep^10 + 
         429011211721788297*ep^9 + 788142926636915931*ep^8 + 
         882399581654092515*ep^7 + 503002260487668978*ep^6 - 5713101262882188*
         ep^5 - 184496065554456216*ep^4 - 87561818231813568*ep^3 - 
         3067248624693120*ep^2 + 6674363573068800*ep + 1242857613312000,
         16777216*ep^11 + 469762048*ep^10 + 5402263552*ep^9 + 32883343360*ep^8
          + 113564975104*ep^7 + 220318400512*ep^6 + 219244658688*ep^5 + 
         84557168640*ep^4) );


id,only intbn*dala^7*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-16*miBN*rat(15116544
         *ep^25 - 161243136*ep^24 - 1506615552*ep^23 + 20027741184*ep^22 + 
         24585729120*ep^21 - 823935442176*ep^20 + 2795513413680*ep^19 + 
         1715940638784*ep^18 - 200612083803975*ep^17 + 1211558672705904*ep^16
          + 9154626232805820*ep^15 - 57834369908763336*ep^14 - 
         448494947086769010*ep^13 + 329990425530430104*ep^12 + 
         11272624044613027380*ep^11 + 45202557254028777864*ep^10 + 
         86081805091705090545*ep^9 + 75627785460694189464*ep^8 - 
         2837918535196828680*ep^7 - 62559786692542255776*ep^6 - 
         39931030542301455024*ep^5 + 5048062290493827840*ep^4 + 
         12307082848011332352*ep^3 + 2240282963258081280*ep^2 - 
         883895358921523200*ep - 237719178805248000,268435456*ep^8 + 
         7516192768*ep^7 + 86436216832*ep^6 + 526133493760*ep^5 + 
         1817039601664*ep^4 + 3525094408192*ep^3 + 3507914539008*ep^2 + 
         1352914698240*ep) + Gam(1,1)^3*M^-16*rat(7558272*ep^22 - 69284160*
         ep^21 - 845896608*ep^20 + 8557958448*ep^19 + 25121929032*ep^18 - 
         360865273068*ep^17 + 825554886354*ep^16 + 23018553911457*ep^15 + 
         409612935180573*ep^14 + 4854736544175711*ep^13 + 12486597245672409*
         ep^12 - 203262008129350875*ep^11 - 2184428082313928061*ep^10 - 
         10602893013763290039*ep^9 - 31252255982156856879*ep^8 - 
         60197233082913291114*ep^7 - 76686020016578500356*ep^6 - 
         62642581670922380424*ep^5 - 29587142078325262656*ep^4 - 
         5218442113063591296*ep^3 + 1614012696910172160*ep^2 + 
         897576105504153600*ep + 118859589402624000,67108864*ep^11 + 
         1879048192*ep^10 + 21609054208*ep^9 + 131533373440*ep^8 + 
         454259900416*ep^7 + 881273602048*ep^6 + 876978634752*ep^5 + 
         338228674560*ep^4) );


id,only intbn*dala^7*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-18*miBN*rat(90699264
         *ep^27 - 136048896*ep^26 - 16471994112*ep^25 + 26247359232*ep^24 + 
         1029784746816*ep^23 - 1953785388384*ep^22 - 11821798414656*ep^21 + 
         31451703633936*ep^20 - 2176923322169946*ep^19 + 6649384311155871*
         ep^18 + 205276275541861839*ep^17 - 578321748894945834*ep^16 - 
         14308894652978231136*ep^15 - 9315029383336648494*ep^14 + 
         690983965966386183714*ep^13 + 5148911488775812662816*ep^12 + 
         18970044518066835858714*ep^11 + 41403081627879826525911*ep^10 + 
         52544833339721954438439*ep^9 + 28269535507104284372466*ep^8 - 
         16783777969883346983400*ep^7 - 34921480093247596422528*ep^6 - 
         15156143725123148129424*ep^5 + 4947848075235662004384*ep^4 + 
         5643031072155927314688*ep^3 + 737938259706696491520*ep^2 - 
         425198770807968460800*ep - 100193994274332672000,2147483648*ep^9 + 
         77309411328*ep^8 + 1172526071808*ep^7 + 9740985827328*ep^6 + 
         48208860413952*ep^5 + 144491289772032*ep^4 + 253669358436352*ep^3 + 
         235329848082432*ep^2 + 86586540687360*ep) + Gam(1,1)^3*M^-18*rat(
         45349632*ep^24 - 8167972608*ep^22 + 372874752*ep^21 + 506192802336*
         ep^20 - 152093427840*ep^19 - 5662880495280*ep^18 + 2712153240864*
         ep^17 - 3829215163823973*ep^16 - 113642091476855640*ep^15 - 
         2138191098711604236*ep^14 - 26679067733202665688*ep^13 - 
         221059293781597756326*ep^12 - 1255437269938050589224*ep^11 - 
         5033740178700919410516*ep^10 - 14511283089558556947624*ep^9 - 
         30253896123152597280333*ep^8 - 45289335604407799858272*ep^7 - 
         47531205020955544324200*ep^6 - 33177141456359875735104*ep^5 - 
         13549882297927013606736*ep^4 - 1801624486034238983424*ep^3 + 
         884343905812130077440*ep^2 + 404637874429788518400*ep + 
         50096997137166336000,536870912*ep^12 + 19327352832*ep^11 + 
         293131517952*ep^10 + 2435246456832*ep^9 + 12052215103488*ep^8 + 
         36122822443008*ep^7 + 63417339609088*ep^6 + 58832462020608*ep^5 + 
         21646635171840*ep^4) );


id,only intbn*dala^8*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-12*miBN*rat(1679616*
         ep^24 - 40310784*ep^23 + 157324032*ep^22 + 3364830720*ep^21 - 
         34128918432*ep^20 + 12009223296*ep^19 + 1516880262384*ep^18 - 
         11019798314208*ep^17 - 2392486627455*ep^16 + 653179013024688*ep^15 - 
         1962455890159044*ep^14 - 29812517272399176*ep^13 + 40872994614441342*
         ep^12 + 1178561385652668120*ep^11 + 5160879809787087828*ep^10 + 
         10701889297573629096*ep^9 + 10404262429732923129*ep^8 + 
         616729535308622232*ep^7 - 8028409335246106920*ep^6 - 
         5696514901754521632*ep^5 + 474146237540206800*ep^4 + 
         1696446119520200448*ep^3 + 335752744805694720*ep^2 - 
         120663702683596800*ep - 33819860938752000,16777216*ep^8 + 469762048*
         ep^7 + 5402263552*ep^6 + 32883343360*ep^5 + 113564975104*ep^4 + 
         220318400512*ep^3 + 219244658688*ep^2 + 84557168640*ep) + Gam(1,1)^3*
         M^-12*rat(839808*ep^21 - 18895680*ep^20 + 51578208*ep^19 + 1722201264
         *ep^18 - 14140510200*ep^17 - 14898381516*ep^16 + 714531348738*ep^15
          - 4387763330487*ep^14 - 92505758659227*ep^13 - 2098194591487401*
         ep^12 - 31410617711700447*ep^11 - 251841528319039155*ep^10 - 
         1200887484686245989*ep^9 - 3660003641008227375*ep^8 - 
         7363678090777500231*ep^7 - 9791587734126974298*ep^6 - 
         8322366352039418628*ep^5 - 4086469336632201864*ep^4 - 
         770085834768679104*ep^3 + 211357616014419840*ep^2 + 
         125153251474406400*ep + 16909930469376000,4194304*ep^11 + 117440512*
         ep^10 + 1350565888*ep^9 + 8220835840*ep^8 + 28391243776*ep^7 + 
         55079600128*ep^6 + 54811164672*ep^5 + 21139292160*ep^4) );


id,only intbn*dala^8*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-14*miBN*rat(5038848*
         ep^24 - 120932352*ep^23 + 471972096*ep^22 + 10094492160*ep^21 - 
         102386755296*ep^20 + 36027669888*ep^19 + 4550640787152*ep^18 - 
         33059394942624*ep^17 - 7177459882365*ep^16 + 1959537039074064*ep^15
          - 5887367670477132*ep^14 - 89437551817197528*ep^13 + 
         122618983843324026*ep^12 + 3535684156958004360*ep^11 + 
         15482639429361263484*ep^10 + 32105667892720887288*ep^9 + 
         31212787289198769387*ep^8 + 1850188605925866696*ep^7 - 
         24085228005738320760*ep^6 - 17089544705263564896*ep^5 + 
         1422438712620620400*ep^4 + 5089338358560601344*ep^3 + 
         1007258234417084160*ep^2 - 361991108050790400*ep - 101459582816256000
         ,67108864*ep^7 + 1744830464*ep^6 + 18119393280*ep^5 + 95294586880*
         ep^4 + 263670726656*ep^3 + 353932148736*ep^2 + 169114337280*ep) + 
         Gam(1,1)^3*M^-14*rat(2519424*ep^21 - 56687040*ep^20 + 154734624*ep^19
          + 5166603792*ep^18 - 42421530600*ep^17 - 44695144548*ep^16 + 
         2143594046214*ep^15 - 13163289991461*ep^14 - 277517275977681*ep^13 - 
         6294583774462203*ep^12 - 94231853135101341*ep^11 - 755524584957117465
         *ep^10 - 3602662454058737967*ep^9 - 10980010923024682125*ep^8 - 
         22091034272332500693*ep^7 - 29374763202380922894*ep^6 - 
         24967099056118255884*ep^5 - 12259408009896605592*ep^4 - 
         2310257504306037312*ep^3 + 634072848043259520*ep^2 + 
         375459754423219200*ep + 50729791408128000,16777216*ep^10 + 436207616*
         ep^9 + 4529848320*ep^8 + 23823646720*ep^7 + 65917681664*ep^6 + 
         88483037184*ep^5 + 42278584320*ep^4) );


id,only intbn*dala^8*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-16*miBN*rat(30233088
         *ep^27 - 367835904*ep^26 - 4630701312*ep^25 + 68409080064*ep^24 + 
         173132157888*ep^23 - 4593825179424*ep^22 + 9440725525056*ep^21 + 
         91464597184176*ep^20 - 1239348796952094*ep^19 + 6203332022517621*
         ep^18 + 76816635798773061*ep^17 - 594845289406647054*ep^16 - 
         5098221203337985824*ep^15 + 16333281541167798726*ep^14 + 
         299831330546954796726*ep^13 + 1419502817349666232416*ep^12 + 
         3435165780989041352286*ep^11 + 4250045595982655488221*ep^10 + 
         1210256954428173317901*ep^9 - 3423421287506181199194*ep^8 - 
         3775645777179955443192*ep^7 - 111761117319556114560*ep^6 + 
         1655155975314530620368*ep^5 + 525195626743365626592*ep^4 - 
         267191151891464017152*ep^3 - 117505771503487879680*ep^2 + 
         14189893467795763200*ep + 7172741757984768000,536870912*ep^9 + 
         19327352832*ep^8 + 293131517952*ep^7 + 2435246456832*ep^6 + 
         12052215103488*ep^5 + 36122822443008*ep^4 + 63417339609088*ep^3 + 
         58832462020608*ep^2 + 21646635171840*ep) + Gam(1,1)^3*M^-16*rat(
         15116544*ep^24 - 161243136*ep^23 - 2534540544*ep^22 + 29994582528*
         ep^21 + 130527508320*ep^20 - 2047190537088*ep^19 + 1523356438896*
         ep^18 + 44879886775008*ep^17 - 1908362238656823*ep^16 - 
         50547606310359432*ep^15 - 886205353450434372*ep^14 - 
         9987945022514554632*ep^13 - 71481221017034524578*ep^12 - 
         337187908169243848440*ep^11 - 1078909282558445893596*ep^10 - 
         2360819161217904535800*ep^9 - 3466706047660611975663*ep^8 - 
         3186815661631723595040*ep^7 - 1419612515176106598264*ep^6 + 
         263145817388399950656*ep^5 + 644238247137196905360*ep^4 + 
         255640354417356247296*ep^3 + 174844747285367040*ep^2 - 
         20842701770035353600*ep - 3586370878992384000,134217728*ep^12 + 
         4831838208*ep^11 + 73282879488*ep^10 + 608811614208*ep^9 + 
         3013053775872*ep^8 + 9030705610752*ep^7 + 15854334902272*ep^6 + 
         14708115505152*ep^5 + 5411658792960*ep^4) );


id,only intbn*dala^8*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-18*miBN*rat(90699264
         *ep^28 - 952342272*ep^27 - 15247554048*ep^26 + 174495306240*ep^25 + 
         793558513728*ep^24 - 11221848109728*ep^23 + 5762270080800*ep^22 + 
         137847889365840*ep^21 - 2459988654875370*ep^20 + 26241694210685385*
         ep^19 + 145431816741459000*ep^18 - 2425808228771702385*ep^17 - 
         9103998912923718630*ep^16 + 119465022493467431730*ep^15 + 
         774819230416416020160*ep^14 - 1069944204921662990610*ep^13 - 
         27370158880915478106630*ep^12 - 129327319034721696202515*ep^11 - 
         320082901311196484294760*ep^10 - 444633964550393305573485*ep^9 - 
         271209597533821906335594*ep^8 + 116132521635702526428072*ep^7 + 
         299137177114105219673328*ep^6 + 141353141601343995169200*ep^5 - 
         38887601604965030724768*ep^4 - 50049341389696649340672*ep^3 - 
         7066643108168236884480*ep^2 + 3726594942997383475200*ep + 
         901745948468994048000,2147483648*ep^9 + 77309411328*ep^8 + 
         1172526071808*ep^7 + 9740985827328*ep^6 + 48208860413952*ep^5 + 
         144491289772032*ep^4 + 253669358436352*ep^3 + 235329848082432*ep^2 + 
         86586540687360*ep) + Gam(1,1)^3*M^-18*rat(45349632*ep^25 - 408146688*
         ep^24 - 8167972608*ep^23 + 73884628224*ep^22 + 502836929568*ep^21 - 
         4707828648864*ep^20 - 4294039644720*ep^19 + 53678077698384*ep^18 - 
         3853624542991749*ep^17 - 79179155002439883*ep^16 - 
         1115412275419903476*ep^15 - 7435347844798227564*ep^14 + 
         19052315817226234866*ep^13 + 734096374096329217710*ep^12 + 
         6265195250741535892500*ep^11 + 30792378518749717747020*ep^10 + 
         100347651682874415248283*ep^9 + 226995729503965575664725*ep^8 + 
         360072815418714654400248*ep^7 + 394603703732240023182696*ep^6 + 
         285044390809311868009200*ep^5 + 120147316195308883477200*ep^4 + 
         17098964280120280928256*ep^3 - 7554457277879382178560*ep^2 - 
         3591643872730930329600*ep - 450872974234497024000,536870912*ep^12 + 
         19327352832*ep^11 + 293131517952*ep^10 + 2435246456832*ep^9 + 
         12052215103488*ep^8 + 36122822443008*ep^7 + 63417339609088*ep^6 + 
         58832462020608*ep^5 + 21646635171840*ep^4) );


id,only intbn*dala^8*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-20*miBN*rat(
         544195584*ep^30 - 140054780160*ep^28 - 4232632320*ep^27 + 
         13253795268480*ep^26 + 1833294145536*ep^25 - 424911131647776*ep^24 - 
         878371617763584*ep^23 - 22122565352089740*ep^22 + 203506148567917440*
         ep^21 + 3496232389190102235*ep^20 - 26487165310177100760*ep^19 - 
         310655735717668257735*ep^18 + 2012289454051748761680*ep^17 + 
         29076072521319304971570*ep^16 - 12861710493694822181520*ep^15 - 
         1899484361547846258112650*ep^14 - 15090210246782293708256640*ep^13 - 
         64248093593662592985309585*ep^12 - 172678845436780702664073240*ep^11
          - 298604947145606411197776579*ep^10 - 304578658611203122644195600*
         ep^9 - 104556716847330730858228140*ep^8 + 145214387713435922164825920
         *ep^7 + 197026718512673533091183280*ep^6 + 65004635072054061924873984
         *ep^5 - 34190198595329347139022144*ep^4 - 28924846690141651178474496*
         ep^3 - 2871133844206328779760640*ep^2 + 2261524130479752787353600*ep
          + 491960527946797400064000,17179869184*ep^10 + 773094113280*ep^9 + 
         14946486190080*ep^8 + 162349763788800*ep^7 + 1087021862879232*ep^6 + 
         4626968267980800*ep^5 + 12432727731077120*ep^4 + 20146832592076800*
         ep^3 + 17636441387433984*ep^2 + 6234230929489920*ep) + Gam(1,1)^3*
         M^-20*rat(272097792*ep^27 + 408146688*ep^26 - 69007023360*ep^25 - 
         108007706880*ep^24 + 6374844379584*ep^23 + 10953535332384*ep^22 - 
         188493102765504*ep^21 - 761913245220912*ep^20 - 12273728886236430*
         ep^19 + 480674018402803107*ep^18 + 23237302549180019229*ep^17 + 
         526562953293857516046*ep^16 + 8122873819146068108448*ep^15 + 
         88763541194931757426890*ep^14 + 696652090474139759220102*ep^13 + 
         3995900348865770815415424*ep^12 + 17013912094570015187070606*ep^11 + 
         54312597425757359425322043*ep^10 + 130365008815103084162924325*ep^9
          + 234078605482539100897704762*ep^8 + 309689066319615921933141192*
         ep^7 + 292916885242242333403363968*ep^6 + 186667653091158042613309008
         *ep^5 + 69525225925124673137974560*ep^4 + 7252737950922177201239808*
         ep^3 - 5051348304416933744494080*ep^2 - 2073686410471238077132800*ep
          - 245980263973398700032000,4294967296*ep^13 + 193273528320*ep^12 + 
         3736621547520*ep^11 + 40587440947200*ep^10 + 271755465719808*ep^9 + 
         1156742066995200*ep^8 + 3108181932769280*ep^7 + 5036708148019200*ep^6
          + 4409110346858496*ep^5 + 1558557732372480*ep^4) );


id,only intbn*dala^9*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-14*miBN*rat(10077696
         *ep^27 - 256981248*ep^26 + 566030592*ep^25 + 35935384320*ep^24 - 
         299599591104*ep^23 - 902371713120*ep^22 + 22232874631104*ep^21 - 
         117254345273712*ep^20 - 207981399667338*ep^19 + 11706826904851671*
         ep^18 - 63134185183443585*ep^17 - 674414842151557386*ep^16 + 
         5110290673106502528*ep^15 + 44199321317598029922*ep^14 - 
         131857061519788187454*ep^13 - 2594417232577918296960*ep^12 - 
         12930538297480930817046*ep^11 - 34092098213381717698737*ep^10 - 
         50619819087367764950313*ep^9 - 34270833211552453869774*ep^8 + 
         10337422367168026942104*ep^7 + 34801473245014912384512*ep^6 + 
         18001576775161556524656*ep^5 - 4063594227768707227488*ep^4 - 
         6103729911565000707840*ep^3 - 937720678490414784000*ep^2 + 
         449704732794960384000*ep + 112426011477319680000,134217728*ep^9 + 
         4831838208*ep^8 + 73282879488*ep^7 + 608811614208*ep^6 + 
         3013053775872*ep^5 + 9030705610752*ep^4 + 15854334902272*ep^3 + 
         14708115505152*ep^2 + 5411658792960*ep) + Gam(1,1)^3*M^-14*rat(
         5038848*ep^24 - 120932352*ep^23 + 109175040*ep^22 + 17894628864*ep^21
          - 121131269856*ep^20 - 617989606272*ep^19 + 9913407771600*ep^18 - 
         43283084425056*ep^17 - 150665554367613*ep^16 + 17834837439120480*
         ep^15 + 446972090124290940*ep^14 + 7863674408189265024*ep^13 + 
         88925385784734434322*ep^12 + 643231342071409422528*ep^11 + 
         3105133650688821099012*ep^10 + 10346876602432971141024*ep^9 + 
         24208678046883771261747*ep^8 + 39791970420537471877728*ep^7 + 
         45119099363804047324440*ep^6 + 33652505597884199194176*ep^5 + 
         14659347350544560209584*ep^4 + 2237597359489033009920*ep^3 - 
         884937532158406944000*ep^2 - 440335555062342912000*ep - 
         56213005738659840000,33554432*ep^12 + 1207959552*ep^11 + 18320719872*
         ep^10 + 152202903552*ep^9 + 753263443968*ep^8 + 2257676402688*ep^7 + 
         3963583725568*ep^6 + 3677028876288*ep^5 + 1352914698240*ep^4) );


id,only intbn*dala^9*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-16*miBN*rat(30233088
         *ep^28 - 700399872*ep^27 - 100776960*ep^26 + 111768367104*ep^25 - 
         647251083072*ep^24 - 4804312277088*ep^23 + 60382021901472*ep^22 - 
         196132913403408*ep^21 - 1444724615917998*ep^20 + 33664610916883647*
         ep^19 - 107454767216369058*ep^18 - 2465183822738777253*ep^17 + 
         10609968124258605882*ep^16 + 168369998664539607462*ep^15 - 
         86175935336178352908*ep^14 - 8706251128372272203058*ep^13 - 
         56952535520488220529858*ep^12 - 192790062722511668815533*ep^11 - 
         390504144755775318742098*ep^10 - 457151233246231716261513*ep^9 - 
         208883565379363096262106*ep^8 + 176766376305220925748264*ep^7 + 
         297615043040589056265552*ep^6 + 113820254742824773990128*ep^5 - 
         46756349329075952715936*ep^4 - 45539271416426249306880*ep^3 - 
         5214930551048022336000*ep^2 + 3485211163996681728000*ep + 
         786982080341237760000,536870912*ep^9 + 19327352832*ep^8 + 
         293131517952*ep^7 + 2435246456832*ep^6 + 12052215103488*ep^5 + 
         36122822443008*ep^4 + 63417339609088*ep^3 + 58832462020608*ep^2 + 
         21646635171840*ep) + Gam(1,1)^3*M^-16*rat(15116544*ep^25 - 327525120*
         ep^24 - 519001344*ep^23 + 54448111872*ep^22 - 238131407520*ep^21 - 
         2701887707808*ep^20 + 25414296070896*ep^19 - 60455398873968*ep^18 - 
         754978254078231*ep^17 + 52449853436788149*ep^16 + 1465760132446716180
         *ep^15 + 26719827855437831652*ep^14 + 321821878211528158134*ep^13 + 
         2552171726707369307838*ep^12 + 13818020346566329254732*ep^11 + 
         52776565362120661116156*ep^10 + 145054170357682111772409*ep^9 + 
         288836657589798814465413*ep^8 + 413901091035174445117416*ep^7 + 
         416791212340280928853608*ep^6 + 279545581236823074987984*ep^5 + 
         109328223532279020496848*ep^4 + 13008368919948010237440*ep^3 - 
         7515569390295877344000*ep^2 - 3250987902652379904000*ep - 
         393491040170618880000,134217728*ep^12 + 4831838208*ep^11 + 
         73282879488*ep^10 + 608811614208*ep^9 + 3013053775872*ep^8 + 
         9030705610752*ep^7 + 15854334902272*ep^6 + 14708115505152*ep^5 + 
         5411658792960*ep^4) );


id,only intbn*dala^9*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-18*miBN*rat(
         181398528*ep^30 - 2176782336*ep^29 - 43419753216*ep^28 + 570236350464
         *ep^27 + 3523516080768*ep^26 - 56125292428800*ep^25 - 43721328214368*
         ep^24 + 1992834170507520*ep^23 - 11927188132929540*ep^22 + 
         97099966469959920*ep^21 + 1066967684793331425*ep^20 - 
         17987624654695426080*ep^19 - 76745786831523706005*ep^18 + 
         1413233702675439882840*ep^17 + 8791327399121406519750*ep^16 - 
         48687343219814557131360*ep^15 - 774843746701108034617470*ep^14 - 
         4142078320855713888888000*ep^13 - 12400221222776467015210755*ep^12 - 
         21942924602624965551873600*ep^11 - 19725860549440923791160393*ep^10
          + 1023747146628972207861816*ep^9 + 20941853598938427514805196*ep^8
          + 16140878980428369740985216*ep^7 - 2330532817521138969771888*ep^6
          - 7930764456686454910861440*ep^5 - 1795826725084616158948032*ep^4 + 
         1338641611389963549987840*ep^3 + 479912017996959037056000*ep^2 - 
         72745033624106425344000*ep - 30790282220370247680000,4294967296*ep^10
          + 193273528320*ep^9 + 3736621547520*ep^8 + 40587440947200*ep^7 + 
         271755465719808*ep^6 + 1156742066995200*ep^5 + 3108181932769280*ep^4
          + 5036708148019200*ep^3 + 4409110346858496*ep^2 + 1558557732372480*
         ep) + Gam(1,1)^3*M^-18*rat(90699264*ep^27 - 952342272*ep^26 - 
         23002341120*ep^25 + 248188458240*ep^24 + 2115999132480*ep^23 - 
         24375987561888*ep^22 - 58455923469120*ep^21 + 865976857887600*ep^20
          - 4452004263580506*ep^19 + 241117197049885401*ep^18 + 
         10478450137738114719*ep^17 + 219393576521890536834*ep^16 + 
         3127002006147879315168*ep^15 + 30736626619840239992670*ep^14 + 
         210483750296547151783938*ep^13 + 1023037280132652966026448*ep^12 + 
         3579004772913173416981722*ep^11 + 9044701741480194809125569*ep^10 + 
         16332012666771647211392247*ep^9 + 20353376684607317583381702*ep^8 + 
         15983471753189181404208792*ep^7 + 5563496361826781294771520*ep^6 - 
         2230572973112762489003664*ep^5 - 3187120132534714527361824*ep^4 - 
         1094489992862816441137920*ep^3 + 34179456828510963648000*ep^2 + 
         95387224401096187392000*ep + 15395141110185123840000,1073741824*ep^13
          + 48318382080*ep^12 + 934155386880*ep^11 + 10146860236800*ep^10 + 
         67938866429952*ep^9 + 289185516748800*ep^8 + 777045483192320*ep^7 + 
         1259177037004800*ep^6 + 1102277586714624*ep^5 + 389639433093120*ep^4)
          );


id,only intbn*dala^9*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-20*miBN*rat(
         544195584*ep^31 - 5441955840*ep^30 - 140054780160*ep^29 + 
         1396315169280*ep^28 + 13296121591680*ep^27 - 130704658539264*ep^26 - 
         443244073103136*ep^25 + 3370739698714176*ep^24 - 13338849174453900*
         ep^23 + 424731802088814840*ep^22 + 1461170903510927835*ep^21 - 
         61449489202078123110*ep^20 - 45784082615897250135*ep^19 + 
         5118846811228431339030*ep^18 + 8953177980801817354770*ep^17 - 
         303622435706887871897220*ep^16 - 1770867256610898036297450*ep^15 + 
         3904633368696168872869860*ep^14 + 86654008874160344097256815*ep^13 + 
         469802090499845227189022610*ep^12 + 1428183507222200615442955821*
         ep^11 + 2681470812844860989333570190*ep^10 + 
         2941229869264700495583727860*ep^9 + 1190781556186743230747107320*ep^8
          - 1255117158621685688557075920*ep^7 - 1905262550054681268986958816*
         ep^6 - 684236549315869966387761984*ep^5 + 312977139263151820211746944
         *ep^4 + 286377333057210183004984320*ep^3 + 30972862572543040584960000
         *ep^2 - 22123280776850730473472000*ep - 4919605279467974000640000,
         17179869184*ep^10 + 773094113280*ep^9 + 14946486190080*ep^8 + 
         162349763788800*ep^7 + 1087021862879232*ep^6 + 4626968267980800*ep^5
          + 12432727731077120*ep^4 + 20146832592076800*ep^3 + 
         17636441387433984*ep^2 + 6234230929489920*ep) + Gam(1,1)^3*M^-20*rat(
         272097792*ep^28 - 2312831232*ep^27 - 73088490240*ep^26 + 582062526720
         *ep^25 + 7454921448384*ep^24 - 52794908463456*ep^23 - 298028456089344
         *ep^22 + 1123017782434128*ep^21 - 4654596434027310*ep^20 + 
         603411307265167407*ep^19 + 18430562365151988159*ep^18 + 
         294189927802057323756*ep^17 + 2857244286207492947988*ep^16 + 
         7534803003471076342410*ep^15 - 190983321475177815048798*ep^14 - 
         2970620555875626776785596*ep^13 - 22945091394087692967083634*ep^12 - 
         115826523519942792445384017*ep^11 - 412760965442470510090296105*ep^10
          - 1069571482668491740731538488*ep^9 - 2031096988505775087043906428*
         ep^8 - 2803973777953916885928047952*ep^7 - 
         2742501199331265291420330672*ep^6 - 1797151304986455752995115520*ep^5
          - 687999521300324554178505792*ep^4 - 77578727813638705756892160*ep^3
          + 48439796633698099367808000*ep^2 + 20490883840738982071296000*ep + 
         2459802639733987000320000,4294967296*ep^13 + 193273528320*ep^12 + 
         3736621547520*ep^11 + 40587440947200*ep^10 + 271755465719808*ep^9 + 
         1156742066995200*ep^8 + 3108181932769280*ep^7 + 5036708148019200*ep^6
          + 4409110346858496*ep^5 + 1558557732372480*ep^4) );


id,only intbn*dala^9*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-22*miBN*rat(
         3265173504*ep^33 + 5986151424*ep^32 - 1139092056576*ep^31 - 
         2242977682176*ep^30 + 152134622102016*ep^29 + 355089430025856*ep^28
          - 8377295133132864*ep^27 - 38635972665742944*ep^26 - 
         111682385417068200*ep^25 + 4749775401491976060*ep^24 + 
         50824575511677228510*ep^23 - 662468457919076761815*ep^22 - 
         4925919138765021735705*ep^21 + 79604237481550312423095*ep^20 + 
         529150600383042646368495*ep^19 - 6876544985586954321801030*ep^18 - 
         74881829624402656771098390*ep^17 + 132829157984336057196468810*ep^16
          + 6632984272334024032107826860*ep^15 + 56093817018756495209903396685
         *ep^14 + 270613847484516190241273219451*ep^13 + 
         863044766465804227583563068531*ep^12 + 
         1880334276482576530334975225031*ep^11 + 
         2724199210222218555416237552256*ep^10 + 
         2291198171149489932796618925964*ep^9 + 370232193063213127794803356224
         *ep^8 - 1433032963863471713850420441456*ep^7 - 
         1486118174588413766104422110976*ep^6 - 366913104033633541162932968640
         *ep^5 + 297541676659176768908501280000*ep^4 + 
         204166847893514411703215232000*ep^3 + 14706173107832430874406400000*
         ep^2 - 16571231247337783311052800000*ep - 
         3373001869828765999104000000,137438953472*ep^11 + 7559142440960*ep^10
          + 181419418583040*ep^9 + 2494517005516800*ep^8 + 21684156006137856*
         ep^7 + 123977495174184960*ep^6 + 469619283287080960*ep^5 + 
         1155792879222784000*ep^4 + 1752838138465615872*ep^3 + 
         1460789158430638080*ep^2 + 498738474359193600*ep) + Gam(1,1)^3*M^-22*
         rat(1632586752*ep^30 + 5441955840*ep^29 - 558934214400*ep^28 - 
         1969685683200*ep^27 + 72322270416000*ep^26 + 288684747429120*ep^25 - 
         3655985720460768*ep^24 - 25072770210817536*ep^23 - 98177725972516644*
         ep^22 + 2216508213736491048*ep^21 - 34986392754763335435*ep^20 - 
         4564854074573730877428*ep^19 - 137853743467274890548651*ep^18 - 
         2604133437300206253876294*ep^17 - 35467394451131609806082022*ep^16 - 
         360051968984974448116739760*ep^15 - 2761725726418491625010723298*
         ep^14 - 16179848001717402246452263116*ep^13 - 
         73038694075893031189824138879*ep^12 - 255425888730744027784705989852*
         ep^11 - 692754312494677954279013646687*ep^10 - 
         1451175895236914272270531408518*ep^9 - 
         2323368639453955158084117536496*ep^8 - 
         2788189985175795082694701090224*ep^7 - 
         2423987247526844498385974231472*ep^6 - 
         1431626115188386827102841642080*ep^5 - 491570270920606085141206896000
         *ep^4 - 38446629150439257838497216000*ep^3 + 
         39165323183982053821728000000*ep^2 + 14750535874174026487142400000*ep
          + 1686500934914382999552000000,34359738368*ep^14 + 1889785610240*
         ep^13 + 45354854645760*ep^12 + 623629251379200*ep^11 + 
         5421039001534464*ep^10 + 30994373793546240*ep^9 + 117404820821770240*
         ep^8 + 288948219805696000*ep^7 + 438209534616403968*ep^6 + 
         365197289607659520*ep^5 + 124684618589798400*ep^4) );


id,only intbn*dala^10*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-16*miBN*rat(
         60466176*ep^30 - 1612431360*ep^29 - 377913600*ep^28 + 346807111680*
         ep^27 - 2263207257216*ep^26 - 21257235025920*ep^25 + 274837242417120*
         ep^24 - 675547557926400*ep^23 - 6650369171439660*ep^22 + 
         166141402971048000*ep^21 - 1249898362101132525*ep^20 - 
         11560157107386058800*ep^19 + 174655330563585433905*ep^18 + 
         758853641445431603400*ep^17 - 13379504199054122134350*ep^16 - 
         85818536286337032638400*ep^15 + 447246901679558774136390*ep^14 + 
         7519334079197023693556400*ep^13 + 41952396153846646193070375*ep^12 + 
         133757584161432798340666800*ep^11 + 264182399292289818744555669*ep^10
          + 307009871367740374400762760*ep^9 + 141613247495193048535500900*
         ep^8 - 117119848838849864095121280*ep^7 - 201344374641185051064765264
         *ep^6 - 78898811098815816397073280*ep^5 + 31099784381727328500742080*
         ep^4 + 31275487346025451988966400*ep^3 + 3706020415997768373120000*
         ep^2 - 2389409004045385405440000*ep - 545672524937663692800000,
         1073741824*ep^10 + 48318382080*ep^9 + 934155386880*ep^8 + 
         10146860236800*ep^7 + 67938866429952*ep^6 + 289185516748800*ep^5 + 
         777045483192320*ep^4 + 1259177037004800*ep^3 + 1102277586714624*ep^2
          + 389639433093120*ep) + Gam(1,1)^3*M^-16*rat(30233088*ep^27 - 
         760866048*ep^26 - 1284906240*ep^25 + 170002333440*ep^24 - 
         868162577472*ep^23 - 11728539034848*ep^22 + 117159826413888*ep^21 - 
         164669307402672*ep^20 - 3340272296540574*ep^19 + 76011722760723435*
         ep^18 - 2491513746435988659*ep^17 - 104425739249840635194*ep^16 - 
         2170146063161190087936*ep^15 - 30954439280527236870822*ep^14 - 
         305858195660392847839386*ep^13 - 2119512352065212886078096*ep^12 - 
         10520557721926696700299458*ep^11 - 38057850379975576655636829*ep^10
          - 101255606360950899277415403*ep^9 - 198069871050662943828013038*
         ep^8 - 281621914643260538576363448*ep^7 - 283275981696297135060813888
         *ep^6 - 190658378109181583087233200*ep^5 - 75119908777301311480153440
         *ep^4 - 9161036236854522654163200*ep^3 + 5114007679690383163200000*
         ep^2 + 2240576841486548113920000*ep + 272836262468831846400000,
         268435456*ep^13 + 12079595520*ep^12 + 233538846720*ep^11 + 
         2536715059200*ep^10 + 16984716607488*ep^9 + 72296379187200*ep^8 + 
         194261370798080*ep^7 + 314794259251200*ep^6 + 275569396678656*ep^5 + 
         97409858273280*ep^4) );


id,only intbn*dala^10*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-18*miBN*rat(
         181398528*ep^31 - 4353564672*ep^30 - 14033191680*ep^29 + 
         1037398026240*ep^28 - 4015164878208*ep^27 - 81877363135488*ep^26 + 
         654453847044000*ep^25 + 172055265557760*ep^24 - 25355487977730180*
         ep^23 + 445221255541626720*ep^22 - 2420563862535013575*ep^21 - 
         44679658218967236600*ep^20 + 431484734831667831315*ep^19 + 
         3673803568844978281440*ep^18 - 34067683465598913575850*ep^17 - 
         364491642451444074990000*ep^16 + 655192414747980061301970*ep^15 + 
         26135977451027541273760320*ep^14 + 186011861095116128127662325*ep^13
          + 736891921715071564566563400*ep^12 + 1862607871168331842959001407*
         ep^11 + 3034488808441539673158733632*ep^10 + 
         2880918713427502140812604780*ep^9 + 781546433444994795998643360*ep^8
          - 1540991914634354065955266032*ep^7 - 1847451430425927857709341952*
         ep^6 - 537891135645344545674360000*ep^5 + 342624737091894983972835840
         *ep^4 + 261321960016196921031091200*ep^3 + 22479936315845990768640000
         *ep^2 - 20752289607176074321920000*ep - 4365380199501309542400000,
         4294967296*ep^10 + 193273528320*ep^9 + 3736621547520*ep^8 + 
         40587440947200*ep^7 + 271755465719808*ep^6 + 1156742066995200*ep^5 + 
         3108181932769280*ep^4 + 5036708148019200*ep^3 + 4409110346858496*ep^2
          + 1558557732372480*ep) + Gam(1,1)^3*M^-18*rat(90699264*ep^28 - 
         2040733440*ep^27 - 9941647104*ep^26 + 499727750400*ep^25 - 
         1244469064896*ep^24 - 42130917724320*ep^23 + 257651166962880*ep^22 + 
         443270689103088*ep^21 - 11338171348843098*ep^20 + 201312989909845713*
         ep^19 - 6866447457222178497*ep^18 - 333209327721009814854*ep^17 - 
         7345844103482295345360*ep^16 - 110224486346871231315954*ep^15 - 
         1165210101225396438484734*ep^14 - 8805402621478781440949376*ep^13 - 
         48517771982301793189523142*ep^12 - 198338012915340303569306151*ep^11
          - 608229622122657311077340841*ep^10 - 1404254464039596025703362338*
         ep^9 - 2429424712335085166353194648*ep^8 - 
         3102803262234975713793349248*ep^7 - 2838182987897921829748210704*ep^6
          - 1750626751205356599138325920*ep^5 - 628442378928974059803717120*
         ep^4 - 57946266855765031743705600*ep^3 + 47633791961982709647360000*
         ep^2 + 18743123519298880450560000*ep + 2182690099750654771200000,
         1073741824*ep^13 + 48318382080*ep^12 + 934155386880*ep^11 + 
         10146860236800*ep^10 + 67938866429952*ep^9 + 289185516748800*ep^8 + 
         777045483192320*ep^7 + 1259177037004800*ep^6 + 1102277586714624*ep^5
          + 389639433093120*ep^4) );


id,only intbn*dala^10*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-20*miBN*rat(
         1088391168*ep^33 - 12516498432*ep^32 - 379697352192*ep^31 + 
         4438323134208*ep^30 + 50824410895872*ep^29 - 608550810910848*ep^28 - 
         2852606574445248*ep^27 + 33795760004639712*ep^26 - 8915434277231160*
         ep^25 + 1040097192158982420*ep^24 + 9172370228045053770*ep^23 - 
         373049190623667620205*ep^22 - 392965575135838051635*ep^21 + 
         44142185150828196973965*ep^20 + 52090040525828646451365*ep^19 - 
         3895958891011281593802210*ep^18 - 19251525717497693634378930*ep^17 + 
         172465495692906065872680270*ep^16 + 2579102469375682983772066020*
         ep^15 + 15345929871288910861275654495*ep^14 + 
         54568824381733432262781986217*ep^13 + 124316137829383998865636171617*
         ep^12 + 174423268559075719611308317677*ep^11 + 
         114574723468384574767561057152*ep^10 - 52349416252403307376710649212*
         ep^9 - 160661226398815533906105475392*ep^8 - 
         92172744320582887161934333392*ep^7 + 30943596466866948798721966848*
         ep^6 + 51812654759336762110112287680*ep^5 + 
         8083259534692832699568326400*ep^4 - 9137615852685725982523008000*ep^3
          - 2767418162684647718561280000*ep^2 + 506987265500079662592000000*ep
          + 187560784581871067136000000,34359738368*ep^11 + 1889785610240*
         ep^10 + 45354854645760*ep^9 + 623629251379200*ep^8 + 5421039001534464
         *ep^7 + 30994373793546240*ep^6 + 117404820821770240*ep^5 + 
         288948219805696000*ep^4 + 438209534616403968*ep^3 + 
         365197289607659520*ep^2 + 124684618589798400*ep) + Gam(1,1)^3*M^-20*
         rat(544195584*ep^30 - 5441955840*ep^29 - 197195316480*ep^28 + 
         1909219507200*ep^27 + 28076020156800*ep^26 - 257786329939200*ep^25 - 
         1798776734030496*ep^24 + 13668477682364928*ep^23 + 16801205434742772*
         ep^22 + 563487234170457816*ep^21 - 26511873396251626665*ep^20 - 
         2174009369866007951916*ep^19 - 58711129194520355658777*ep^18 - 
         1022688307479483849573378*ep^17 - 12785324706791763423955074*ep^16 - 
         116933715223245694662254640*ep^15 - 790429808256099743609612406*ep^14
          - 3989248081829961795885698532*ep^13 - 15138833026342484473508836773
         *ep^12 - 43247930893130717673725029284*ep^11 - 
         92326607175740735649417894549*ep^10 - 144343650007219624240596034626*
         ep^9 - 158057715471007245858249519312*ep^8 - 
         108536955622053450212878758288*ep^7 - 28738061220971935003366817424*
         ep^6 + 20378296834934579875843745760*ep^5 + 
         21678548114914884887383036800*ep^4 + 6580389789925286272530624000*
         ep^3 - 408594943411848079407360000*ep^2 - 612985136531959376640000000
         *ep - 93780392290935533568000000,8589934592*ep^14 + 472446402560*
         ep^13 + 11338713661440*ep^12 + 155907312844800*ep^11 + 
         1355259750383616*ep^10 + 7748593448386560*ep^9 + 29351205205442560*
         ep^8 + 72237054951424000*ep^7 + 109552383654100992*ep^6 + 
         91299322401914880*ep^5 + 31171154647449600*ep^4) );


id,only intbn*dala^10*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-22*miBN*rat(
         3265173504*ep^34 - 29930757120*ep^33 - 1204939722240*ep^32 + 
         10287034940160*ep^31 + 176807376605952*ep^30 - 1318391413096320*ep^29
          - 12283278863417280*ep^28 + 53514273798718560*ep^27 + 
         313313313906104184*ep^26 + 5978281641079726260*ep^25 - 
         1422953904734508150*ep^24 - 1221538788547526275425*ep^23 + 
         2361233898344822644260*ep^22 + 133789348007965551515850*ep^21 - 
         346496011914010790285550*ep^20 - 12697201589800423431854475*ep^19 + 
         760165217053840768712940*ep^18 + 956529283852765281678551100*ep^17 + 
         5171863534506327402946669950*ep^16 - 16869009976917769143282698775*
         ep^15 - 346418139721805257067664144084*ep^14 - 
         2113707555863873865070442345430*ep^13 - 
         7613158154641269973084218528810*ep^12 - 
         17959477831086123278268489923085*ep^11 - 
         27674993141294914176781994148852*ep^10 - 
         24832947689581176132968004829380*ep^9 - 
         5505587087558816119593257359920*ep^8 + 
         14277244427909775086250202745040*ep^7 + 
         15980386816438917885985710252096*ep^6 + 
         4333585821029145721700763935040*ep^5 - 
         3068791595357430046290298848000*ep^4 - 
         2231129153720826097860961152000*ep^3 - 178339135433494522929523200000
         *ep^2 + 178910541850886850422476800000*ep + 
         37103020568116425990144000000,137438953472*ep^11 + 7559142440960*
         ep^10 + 181419418583040*ep^9 + 2494517005516800*ep^8 + 
         21684156006137856*ep^7 + 123977495174184960*ep^6 + 469619283287080960
         *ep^5 + 1155792879222784000*ep^4 + 1752838138465615872*ep^3 + 
         1460789158430638080*ep^2 + 498738474359193600*ep) + Gam(1,1)^3*M^-22*
         rat(1632586752*ep^31 - 12516498432*ep^30 - 618795728640*ep^29 + 
         4178590675200*ep^28 + 93988812931200*ep^27 - 506860227146880*ep^26 - 
         6831517942181088*ep^25 + 15143072714250912*ep^24 + 177622746346476252
         *ep^23 + 3296463199434174132*ep^22 - 59367983105864736963*ep^21 - 
         4180003754271334187643*ep^20 - 87640348646963850896943*ep^19 - 
         1087742259160182457841133*ep^18 - 6821926640829341013442788*ep^17 + 
         30089369977473259750162482*ep^16 + 1198845932416227304273414062*ep^15
          + 14199134988886005628665693162*ep^14 + 
         104939633942998393521150755397*ep^13 + 547999746104079315303359537817
         *ep^12 + 2116930463543506351352752241685*ep^11 + 
         6169121542204543224798618705039*ep^10 + 
         13639566208152101836891727957202*ep^9 + 
         22768865048817711656230591811232*ep^8 + 
         28246102589406901411255737760992*ep^7 + 
         25232233607606902655142874904112*ep^6 + 
         15256316996151649012990051166880*ep^5 + 
         5368826350976227678714778640000*ep^4 + 462078243838813890045197376000
         *ep^3 - 416068019149628565551865600000*ep^2 - 
         160569393680999908359014400000*ep - 18551510284058212995072000000,
         34359738368*ep^14 + 1889785610240*ep^13 + 45354854645760*ep^12 + 
         623629251379200*ep^11 + 5421039001534464*ep^10 + 30994373793546240*
         ep^9 + 117404820821770240*ep^8 + 288948219805696000*ep^7 + 
         438209534616403968*ep^6 + 365197289607659520*ep^5 + 
         124684618589798400*ep^4) );


id,only intbn*dala^10*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-24*miBN*rat(
         19591041024*ep^36 + 78364164096*ep^35 - 8941133445120*ep^34 - 
         38121989050368*ep^33 + 1602649207551744*ep^32 + 7817986923816960*
         ep^31 - 129467506239422208*ep^30 - 963239522958484992*ep^29 + 
         1367431245248829552*ep^28 + 103073015495689614144*ep^27 + 
         680361892509151089288*ep^26 - 13475107885183654280400*ep^25 - 
         68677506788506593048135*ep^24 + 1941102306018155084843520*ep^23 + 
         4722868610817737013500130*ep^22 - 256464716144965835096428740*ep^21
          - 917392642218552511646333325*ep^20 + 26561732989954897972211394420*
         ep^19 + 238982173488701029493059143960*ep^18 - 
         828001361557521028387706989320*ep^17 - 
         28709516473350672436395672527649*ep^16 - 
         257846595976727315185371223979736*ep^15 - 
         1388779447897560847583604992826150*ep^14 - 
         5123784539350667201184159587985492*ep^13 - 
         13451397611907859164197830908645579*ep^12 - 
         25026255065374177655693288600413020*ep^11 - 
         31351059224744618890074255595027452*ep^10 - 
         21937157132011915543113003030431568*ep^9 + 
         629409471015685201612703222680368*ep^8 + 
         17265598682329160347629266372903616*ep^7 + 
         14676867679660186526492140683087552*ep^6 + 
         2573238779464291113104111553050880*ep^5 - 
         3297628271940500292332122482048000*ep^4 - 
         1926207101275417321994823289344000*ep^3 - 
         91703550183707227240503705600000*ep^2 + 
         162347768650437671194487193600000*ep + 
         31228469227906297060589568000000,1099511627776*ep^12 + 72567767433216
         *ep^11 + 2116559883468800*ep^10 + 35921044879441920*ep^9 + 
         392990744534581248*ep^8 + 2900025689933611008*ep^7 + 
         14666973841624924160*ep^6 + 50572839963045396480*ep^5 + 
         115732478479329918976*ep^4 + 165936069452419301376*ep^3 + 
         132539353736769699840*ep^2 + 43888985743609036800*ep) + Gam(1,1)^3*
         M^-24*rat(9795520512*ep^33 + 53875362816*ep^32 - 4375060397568*ep^31
          - 25650522802944*ep^30 + 756340104337920*ep^29 + 5051331531768960*
         ep^28 - 56002249062767808*ep^27 - 566046056811327840*ep^26 - 
         259825949649945720*ep^25 + 50888665987139518236*ep^24 + 
         419635528008516070362*ep^23 + 5178537392413495846833*ep^22 + 
         878432481614586229450569*ep^21 + 36576770887604900460485385*ep^20 + 
         865654877912942089764853821*ep^19 + 14426228026554761017690799832*
         ep^18 + 181579749950473160176206018018*ep^17 + 
         1769679535544441920038295855134*ep^16 + 
         13507323970046566289374545162660*ep^15 + 
         81336709192397795743504933306161*ep^14 + 
         388422719969888862949987000910925*ep^13 + 
         1475245603492303444675240904873037*ep^12 + 
         4456127751847278752594021688673317*ep^11 + 
         10667458203783993867270965849789814*ep^10 + 
         20080232645001010983344576102398416*ep^9 + 
         29320637164020458822618505111258432*ep^8 + 
         32483599122039644241107933778531216*ep^7 + 
         26318130249039674084157436311186144*ep^6 + 
         14562136868081205363021190227807360*ep^5 + 
         4645881605981242416226161312576000*ep^4 + 
         249684095692602874822867875072000*ep^3 - 
         401939334360241972475927961600000*ep^2 - 
         141028450345372571630040268800000*ep - 
         15614234613953148530294784000000,274877906944*ep^15 + 18141941858304*
         ep^14 + 529139970867200*ep^13 + 8980261219860480*ep^12 + 
         98247686133645312*ep^11 + 725006422483402752*ep^10 + 
         3666743460406231040*ep^9 + 12643209990761349120*ep^8 + 
         28933119619832479744*ep^7 + 41484017363104825344*ep^6 + 
         33134838434192424960*ep^5 + 10972246435902259200*ep^4) );


id,only intbn*dala^11*x3^1*x4^1*x5^1*x6^1 =  + int0 * ( M^-18*miBN*rat(
         362797056*ep^33 - 9976919040*ep^32 - 35140925952*ep^31 + 
         3106868016384*ep^30 - 13856573339136*ep^29 - 331022949713280*ep^28 + 
         2857385683703232*ep^27 + 5752804708338336*ep^26 - 132182086669090920*
         ep^25 + 1755344907442734300*ep^24 - 18789247456996290930*ep^23 - 
         163261528249375164015*ep^22 + 4025834205111409209255*ep^21 + 
         5710593787674411905775*ep^20 - 461638380023767243273185*ep^19 - 
         597446399360175994606230*ep^18 + 40853137355505481288323690*ep^17 + 
         209650190205785897008618650*ep^16 - 1774167676651744318376211780*
         ep^15 - 27725136672381183891963140715*ep^14 - 
         170982718619873227930357646661*ep^13 - 638290617192817152228045197685
         *ep^12 - 1567903964537249006522383939113*ep^11 - 
         2521985656123055253776107646304*ep^10 - 
         2387972624350100591131608392244*ep^9 - 657627764194857865812368454720
         *ep^8 + 1270210935604127041721345635728*ep^7 + 
         1542304566983518518799778252544*ep^6 + 459049491577079572197680688960
         *ep^5 - 283426405806897111893248224000*ep^4 - 
         220447967361280413759128448000*ep^3 - 19603431363406186313740800000*
         ep^2 + 17475378519951316313395200000*ep + 
         3705756097386556071936000000,8589934592*ep^11 + 472446402560*ep^10 + 
         11338713661440*ep^9 + 155907312844800*ep^8 + 1355259750383616*ep^7 + 
         7748593448386560*ep^6 + 29351205205442560*ep^5 + 72237054951424000*
         ep^4 + 109552383654100992*ep^3 + 91299322401914880*ep^2 + 
         31171154647449600*ep) + Gam(1,1)^3*M^-18*rat(181398528*ep^30 - 
         4716361728*ep^29 - 24372907776*ep^28 + 1507804720128*ep^27 - 
         4639286669184*ep^26 - 170351463460608*ep^25 + 1151641765269792*ep^24
          + 4472292047567616*ep^23 - 56533282199472132*ep^22 + 
         780754720731217416*ep^21 - 8213591161117198899*ep^20 + 
         257391591278589679860*ep^19 + 23648524272629904113397*ep^18 + 
         639673407078507739153530*ep^17 + 11148476086212431057851050*ep^16 + 
         139775669927962679344371264*ep^15 + 1286799608013892190212796766*
         ep^14 + 8801676948514298126904722484*ep^13 + 
         45262389698551305081246632409*ep^12 + 176636920376742898129922213388*
         ep^11 + 525603337455638697652262898657*ep^10 + 
         1191303649544349635799267525786*ep^9 + 
         2040273931113120697075579985424*ep^8 + 
         2594875963154625733380373281744*ep^7 + 
         2373478051966094228365562341968*ep^6 + 
         1468311923553874713389487561120*ep^5 + 530195550412751741731007376000
         *ep^4 + 50111551794437910947663424000*ep^3 - 
         39905442740936445925190400000*ep^2 - 15840388446633223961241600000*ep
          - 1852878048693278035968000000,2147483648*ep^14 + 118111600640*ep^13
          + 2834678415360*ep^12 + 38976828211200*ep^11 + 338814937595904*ep^10
          + 1937148362096640*ep^9 + 7337801301360640*ep^8 + 18059263737856000*
         ep^7 + 27388095913525248*ep^6 + 22824830600478720*ep^5 + 
         7792788661862400*ep^4) );


id,only intbn*dala^11*x3^1*x4^1*x5^1*x6^2 =  + int0 * ( M^-20*miBN*rat(
         1088391168*ep^33 - 29930757120*ep^32 - 105422777856*ep^31 + 
         9320604049152*ep^30 - 41569720017408*ep^29 - 993068849139840*ep^28 + 
         8572157051109696*ep^27 + 17258414125015008*ep^26 - 396546260007272760
         *ep^25 + 5266034722328202900*ep^24 - 56367742370988872790*ep^23 - 
         489784584748125492045*ep^22 + 12077502615334227627765*ep^21 + 
         17131781363023235717325*ep^20 - 1384915140071301729819555*ep^19 - 
         1792339198080527983818690*ep^18 + 122559412066516443864971070*ep^17
          + 628950570617357691025855950*ep^16 - 5322503029955232955128635340*
         ep^15 - 83175410017143551675889422145*ep^14 - 
         512948155859619683791072939983*ep^13 - 
         1914871851578451456684135593055*ep^12 - 
         4703711893611747019567151817339*ep^11 - 
         7565956968369165761328322938912*ep^10 - 
         7163917873050301773394825176732*ep^9 - 
         1972883292584573597437105364160*ep^8 + 
         3810632806812381125164036907184*ep^7 + 
         4626913700950555556399334757632*ep^6 + 
         1377148474731238716593042066880*ep^5 - 850279217420691335679744672000
         *ep^4 - 661343902083841241277385344000*ep^3 - 
         58810294090218558941222400000*ep^2 + 52426135559853948940185600000*ep
          + 11117268292159668215808000000,34359738368*ep^10 + 1786706395136*
         ep^9 + 39994735460352*ep^8 + 503645044998144*ep^7 + 3910103866540032*
         ep^6 + 19264062193926144*ep^5 + 59612634239991808*ep^4 + 
         110110317085720576*ep^3 + 107878583359242240*ep^2 + 41561539529932800
         *ep) + Gam(1,1)^3*M^-20*rat(544195584*ep^30 - 14149085184*ep^29 - 
         73118723328*ep^28 + 4523414160384*ep^27 - 13917860007552*ep^26 - 
         511054390381824*ep^25 + 3454925295809376*ep^24 + 13416876142702848*
         ep^23 - 169599846598416396*ep^22 + 2342264162193652248*ep^21 - 
         24640773483351596697*ep^20 + 772174773835769039580*ep^19 + 
         70945572817889712340191*ep^18 + 1919020221235523217460590*ep^17 + 
         33445428258637293173553150*ep^16 + 419327009783888038033113792*ep^15
          + 3860398824041676570638390298*ep^14 + 26405030845542894380714167452
         *ep^13 + 135787169095653915243739897227*ep^12 + 
         529910761130228694389766640164*ep^11 + 
         1576810012366916092956788695971*ep^10 + 
         3573910948633048907397802577358*ep^9 + 
         6120821793339362091226739956272*ep^8 + 
         7784627889463877200141119845232*ep^7 + 
         7120434155898282685096687025904*ep^6 + 
         4404935770661624140168462683360*ep^5 + 
         1590586651238255225193022128000*ep^4 + 150334655383313732842990272000
         *ep^3 - 119716328222809337775571200000*ep^2 - 
         47521165339899671883724800000*ep - 5558634146079834107904000000,
         8589934592*ep^13 + 446676598784*ep^12 + 9998683865088*ep^11 + 
         125911261249536*ep^10 + 977525966635008*ep^9 + 4816015548481536*ep^8
          + 14903158559997952*ep^7 + 27527579271430144*ep^6 + 
         26969645839810560*ep^5 + 10390384882483200*ep^4) );


id,only intbn*dala^11*x3^1*x4^1*x5^2*x6^2 =  + int0 * ( M^-22*miBN*rat(
         6530347008*ep^36 - 69657034752*ep^35 - 3155971590144*ep^34 + 
         32295226466304*ep^33 + 621257261315328*ep^32 - 5899637346527232*ep^31
          - 61532867480502528*ep^30 + 464630634966491136*ep^29 + 
         2849463030249841104*ep^28 + 7899875611597626624*ep^27 - 
         49671393255657890568*ep^26 - 6378954856710730860960*ep^25 + 
         13398441994486241026155*ep^24 + 971921090092605375586680*ep^23 - 
         3112089475827006613179450*ep^22 - 118493718656081179806214380*ep^21
          + 178456412405549065117757625*ep^20 + 12928576529926989955533262500*
         ep^19 + 51357252193183188947077507080*ep^18 - 
         736348566454518210837283803240*ep^17 - 
         10743802498928011225946390515683*ep^16 - 
         70344199341184017344239708962048*ep^15 - 
         289146807910277486205302391651666*ep^14 - 
         802585902555335180721520055691324*ep^13 - 
         1500475489693176993921560086508673*ep^12 - 
         1735360376149619291590858289922588*ep^11 - 
         792122155762666155749366944651092*ep^10 + 
         861875965069279590708077430693744*ep^9 + 
         1567473420705519665929652835173136*ep^8 + 
         679650678755291096724109204506816*ep^7 - 
         391167436967750270217231355741632*ep^6 - 
         450539544966822700476227212481280*ep^5 - 
         43566214035691770303782871936000*ep^4 + 
         83073254481496619584128382464000*ep^3 + 
         21728201744345809372698931200000*ep^2 - 
         4700900557071681833867673600000*ep - 1562418227343467765956608000000,
         274877906944*ep^12 + 18141941858304*ep^11 + 529139970867200*ep^10 + 
         8980261219860480*ep^9 + 98247686133645312*ep^8 + 725006422483402752*
         ep^7 + 3666743460406231040*ep^6 + 12643209990761349120*ep^5 + 
         28933119619832479744*ep^4 + 41484017363104825344*ep^3 + 
         33134838434192424960*ep^2 + 10972246435902259200*ep) + Gam(1,1)^3*
         M^-22*rat(3265173504*ep^33 - 29930757120*ep^32 - 1617984170496*ep^31
          + 13639823933184*ep^30 + 329206130081280*ep^29 - 2421519926590080*
         ep^28 - 34126573330414656*ep^27 + 175256073992929248*ep^26 + 
         1673497673175755736*ep^25 + 6866456640418241748*ep^24 - 
         14130684001693421730*ep^23 + 2403981176225770713267*ep^22 + 
         434580339878980234447779*ep^21 + 16111788400890790277684931*ep^20 + 
         348737351403099717907704783*ep^19 + 5359830692604933215354660640*
         ep^18 + 61845170126926105139734779462*ep^17 + 
         545018194451206600608026045034*ep^16 + 
         3699105742815229921171766450796*ep^15 + 
         19452419578724094555031758616803*ep^14 + 
         79539058885302479424939143452911*ep^13 + 
         252830656157242036180848891665967*ep^12 + 
         621352122474388266900759690892935*ep^11 + 
         1165155574870870354673143627692762*ep^10 + 
         1624490810404605926495030453479488*ep^9 + 
         1598239589632368871732262438339520*ep^8 + 
         973645086008118747119935544776368*ep^7 + 
         183590740404949066945336300440096*ep^6 - 
         226144901026499639093359940684160*ep^5 - 
         197268395740520533031240593344000*ep^4 - 
         53664552886067218990773004032000*ep^3 + 
         4981538273508823637957068800000*ep^2 + 
         5345085214277487468350668800000*ep + 781209113671733882978304000000,
         68719476736*ep^15 + 4535485464576*ep^14 + 132284992716800*ep^13 + 
         2245065304965120*ep^12 + 24561921533411328*ep^11 + 181251605620850688
         *ep^10 + 916685865101557760*ep^9 + 3160802497690337280*ep^8 + 
         7233279904958119936*ep^7 + 10371004340776206336*ep^6 + 
         8283709608548106240*ep^5 + 2743061608975564800*ep^4) );


id,only intbn*dala^11*x3^1*x4^2*x5^2*x6^2 =  + int0 * ( M^-24*miBN*rat(
         19591041024*ep^37 - 156728328192*ep^36 - 9881503414272*ep^35 + 
         69171612291072*ep^34 + 2060113076156160*ep^33 - 11413803566803968*
         ep^32 - 223283349325225728*ep^31 + 590370551914581504*ep^30 + 
         12926305520750649456*ep^29 + 86663840552703659520*ep^28 - 
         556514293439124280440*ep^27 - 21639450595293467351856*ep^26 + 
         93023787833697258316665*ep^25 + 2765232387480234201421140*ep^24 - 
         18570359061400124004622110*ep^23 - 313139139474778679258430300*ep^22
          + 2160183951521037509510811555*ep^21 + 37570444696577528111967394320
         *ep^20 - 79758622390757746173477589080*ep^19 - 
         3695787443421933382304416716840*ep^18 - 
         18773500134660420095743188655809*ep^17 + 
         86667601703480754051376846352052*ep^16 + 
         1705379703823166934640849694930682*ep^15 + 
         11541568835420062969819100325928308*ep^14 + 
         48034016860300147250012084147180325*ep^13 + 1363905162775201323146806\
         82303333928*ep^12 + 268964001559745512978245207609928788*ep^11 + 
         354275553564923511137778064109897856*ep^10 + 263875295055158671718968\
         739587859184*ep^9 + 9712685030140937928276827700739200*ep^8 - 
         192510316508289737645059055791755840*ep^7 - 1735491733764579472048015\
         76643999744*ep^6 - 34176493625511993649581461118658560*ep^5 + 
         37645332162010586185990646495232000*ep^4 + 
         23022781665121300636697375766528000*ep^3 + 
         1262790370854924398080531660800000*ep^2 - 
         1916944754577345757273256755200000*ep - 
         374741630734875564727074816000000,1099511627776*ep^12 + 
         72567767433216*ep^11 + 2116559883468800*ep^10 + 35921044879441920*
         ep^9 + 392990744534581248*ep^8 + 2900025689933611008*ep^7 + 
         14666973841624924160*ep^6 + 50572839963045396480*ep^5 + 
         115732478479329918976*ep^4 + 165936069452419301376*ep^3 + 
         132539353736769699840*ep^2 + 43888985743609036800*ep) + Gam(1,1)^3*
         M^-24*rat(9795520512*ep^34 - 63670883328*ep^33 - 5021564751360*ep^32
          + 26850201967872*ep^31 + 1064146377973248*ep^30 - 4024749720286080*
         ep^29 - 116618227443995328*ep^28 + 105980931941885856*ep^27 + 
         6532726732085988360*ep^26 + 54006577382938866876*ep^25 - 
         191028463837158148470*ep^24 + 142911056311303002489*ep^23 + 
         816290032905624279288573*ep^22 + 26035581108229865707078557*ep^21 + 
         426733627261683284239029201*ep^20 + 4038369491599455940512553980*
         ep^19 + 8465013631816027963916420034*ep^18 - 
         409277463861236002076176361082*ep^17 - 
         7728830456486736751085005098948*ep^16 - 
         80751178448160999728989608645759*ep^15 - 
         587617790338884685972072198763007*ep^14 - 
         3185827036146362910724603106058063*ep^13 - 
         13246819490060362583508869169803127*ep^12 - 
         42806074818383351163857294414289990*ep^11 - 1079292658004069154239070\
         14095079352*ep^10 - 211642154575991672977516408117522560*ep^9 - 
         319364046846205861630314127556569968*ep^8 - 3634850592154360568091377\
         69031188448*ep^7 - 301255426120394883646868045506426368*ep^6 - 
         170099760810993221940028121421112320*ep^5 - 
         55500895176082306119891067875840000*ep^4 - 
         3398148482671476470350342462464000*ep^3 + 
         4682243561977531098081095270400000*ep^2 + 
         1676727169530517711030188441600000*ep + 
         187370815367437782363537408000000,274877906944*ep^15 + 18141941858304
         *ep^14 + 529139970867200*ep^13 + 8980261219860480*ep^12 + 
         98247686133645312*ep^11 + 725006422483402752*ep^10 + 
         3666743460406231040*ep^9 + 12643209990761349120*ep^8 + 
         28933119619832479744*ep^7 + 41484017363104825344*ep^6 + 
         33134838434192424960*ep^5 + 10972246435902259200*ep^4) );


id,only intbn*dala^11*x3^2*x4^2*x5^2*x6^2 =  + int0 * ( M^-26*miBN*rat(
         117546246144*ep^39 + 764050599936*ep^38 - 68163762069504*ep^37 - 
         473930138585088*ep^36 + 15772167650824704*ep^35 + 124707552876382464*
         ep^34 - 1717847922431953152*ep^33 - 18955118424913761024*ep^32 + 
         49485051068399701920*ep^31 + 2131544589272534194704*ep^30 + 
         9185804495047780880256*ep^29 - 255521678386522145190408*ep^28 - 
         1033101614468459222211618*ep^27 + 38212718159020864066545591*ep^26 - 
         16597308281441078659345245*ep^25 - 6072240237232074709835032020*ep^24
          + 11834329595021542071248282820*ep^23 + 
         936009633646094203530579046725*ep^22 + 918406534159715755057076201985
         *ep^21 - 119395078208054202227541840073530*ep^20 - 
         924673004277439424147840351921454*ep^19 + 
         5244584967371225707121888213422329*ep^18 + 15094091487726711436382903\
         9697576309*ep^17 + 1436020004178076261229808028632643248*ep^16 + 
         8535977551588009138511855199394163256*ep^15 + 35762163649750460888461\
         249097983627491*ep^14 + 109874861296835618027897543713240626327*ep^13
          + 249146219179515273759925413361790407974*ep^12 + 408201825433476028\
         733345559718288623060*ep^11 + 450881398248611050810239394921454268216
         *ep^10 + 261506715636871759774169268106156833264*ep^9 - 6530846200873\
         5668875412702439676141152*ep^8 - 257112687634003617856991165159301920\
         832*ep^7 - 186063735309283830597840339973817187456*ep^6 - 20461025267\
         393573399921080676517450240*ep^5 + 4624560640156279339067576364652723\
         2000*ep^4 + 23690794181716538145445788271475712000*ep^3 + 59415359615\
         1716364731110464921600000*ep^2 - 2073992136723656453875266817228800000
         *ep - 379877015415849782690317860864000000,8796093022208*ep^13 + 
         686095255732224*ep^12 + 23898984741339136*ep^11 + 490558107848540160*
         ep^10 + 6592346264703074304*ep^9 + 60927316994788687872*ep^8 + 
         395738256966626050048*ep^7 + 1812612208500355891200*ep^6 + 
         5780852464286997413888*ep^5 + 12437806489635026632704*ep^4 + 
         16990177497326410530816*ep^3 + 13074889844678763479040*ep^2 + 
         4213342631386467532800*ep) + Gam(1,1)^3*M^-26*rat(58773123072*ep^36
          + 470184984576*ep^35 - 33288443873280*ep^34 - 286838961979392*ep^33
          + 7404599707925760*ep^32 + 73395356204132352*ep^31 - 
         736686933144028416*ep^30 - 10556095513519319040*ep^29 + 
         7477919438265561168*ep^28 + 1070267874588428909184*ep^27 + 
         6275512254014221315320*ep^26 - 117177451528637290295232*ep^25 - 
         2847806458414012074010077*ep^24 - 194660375948344105847074152*ep^23
          - 9895504639437722640727314702*ep^22 - 
         293100939355154997812128084512*ep^21 - 
         5966875760161228205460501912867*ep^20 - 
         91277399940832290749003573387832*ep^19 - 
         1092863670607994352708051499377672*ep^18 - 
         10431330020087185059014145274498512*ep^17 - 
         80103988363760414000723990259535779*ep^16 - 4975369127359427462402128\
         64851406808*ep^15 - 2507612103080672492568368537076876102*ep^14 - 
         10269952146214799100609751372828650240*ep^13 - 3415565109366136728063\
         6220547943673437*ep^12 - 91956963747734557809497699874822574344*ep^11
          - 199184801601537898976504689140428029116*ep^10 - 343631779683690525\
         645301608721126419792*ep^9 - 464881702323681656012599674034987545840*
         ep^8 - 481507482312589606890150346083647006976*ep^7 - 367268995906878\
         886098487666223845396032*ep^6 - 1919353066676336628222634051939928332\
         80*ep^5 - 57229067783694783072253262934506112000*ep^4 - 1734097151163\
         089072756617456928256000*ep^3 + 5340036739755492820165513475174400000
         *ep^2 + 1765093681242206977094075975270400000*ep + 189938507707924891\
         345158930432000000,2199023255552*ep^16 + 171523813933056*ep^15 + 
         5974746185334784*ep^14 + 122639526962135040*ep^13 + 
         1648086566175768576*ep^12 + 15231829248697171968*ep^11 + 
         98934564241656512512*ep^10 + 453153052125088972800*ep^9 + 
         1445213116071749353472*ep^8 + 3109451622408756658176*ep^7 + 
         4247544374331602632704*ep^6 + 3268722461169690869760*ep^5 + 
         1053335657846616883200*ep^4) );
#endprocedure

#procedure reduceBNBN
*
* reduce BN to BM and simpler integrals
*

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

        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );
        multiply replace_(x3,x4,x4,x3);
        endif;
        endif;

        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );
        multiply replace_(x3,x5,x5,x3);
        endif;
        endif;
        
* #include expandnomdeno

        #call redBNn3456(p1,p2,x3,x4,x5,x6)
        .sort
        
        if ( count(intbn,1) );        
        if ( (count(p1.p1,1)==0) && (count(p2.p2,1)==0) );
        multiply replace_(x3,x6,x6,x3);
        endif;
        endif;

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


* Check if table is too small:
        
        if ( count(intbn,1) );   
        
        if ( count(dala,1) > 11) exit "ERROR: Table is too small for dala > 11 in BN reduction"; 
*         if (count(dala,1) > 11) multiply 1/(1-1); 

        #call BNdExact
        
* Replace temporarily x6 by s6m, in order not to touch
* terms, which are not listed in the previous table.

        if ( (count(dala,1)!=0) ) id x6=s6m;
        endif;        

        .sort

        if ( count(intbn,1) );                
        #call BNd0Exact

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
        #ifdef `REDBNTAB'
                #call reduceBNBN
                #else                
                #call reduceBNnotab
        #endif

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




************************************************************
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
*         .sort:Before match;

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

        #call Conv2exact        
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





************************************************************
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
        
        #call Conv2exact
        
        
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





************************************************************
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

*         Print+s;
*         .end

        
#endprocedure        


#procedure symBM2 (p1,p2,p3,p4,p5,x6)

        if( count(intbm2,1));        
        if (count(`p3'.`p3',1) >= count(`p1'.`p1',1))
        multiply replace_(`p1',`p3',`p3',`p1',`p4',`p5',`p5',`p4');
        endif;
#endprocedure




************************************************************
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
*         #call triown(p2,p3,p5,p1,p4)
        #call ntriangle(p2,p3,p5,p1,p4)        
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



#procedure nomgm3(v1,v2,v3,x,y,z,in)
* Apply only to topo "in"
if(count(int`in',1)) id  `v3' = -`v1'-`v2';

#call ACCU{nomgm3 1}

if(count(int`in',1)) id  `v1'.`v1' = 1/`x' - M^2;

#call ACCU{nomgm3 2}

if(count(int`in',1)) id  `v2'.`v2' = 1/`y' - M^2;

#call ACCU{nomgm3 3}

if(count(int`in',1)) id  `v1'.`v2' = (1/(`z') - 1/`x' - 1/`y' +2*M^2)/2;

#call ACCU{nomgm3 4}

#endprocedure



************************************************************
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
        #call nomgm3(p5,p6,p3,x5,x6,1/p3.p3,MM0)
        
        .sort

*         #call two110(x5,s5m1,x6,s6m1,1/p3.p3,e3)
        #call TadpoleMM0(x5,x6,p3,MM0,0)
        .sort

        #call ACCU(gm3)
* #include expandnomdeno
        
#endprocedure        


************************************************************
#procedure topm2
*
* topm2
* (simple topology)
*
        
        #message this is topm2
*         
* Procedure m2 modified        
* IntOne can deal only with numerator containing first momentum        
*         
        if( count(intm2,1));        
***         id p1=-p1;
        id p1=p2-p3;
        endif;        
        #call IntOne(p2,p1,p3,m2,M00)        
        .sort

        if( count(intM00,1));        
        id p5=p6-p3;
        endif;        
        
        #call IntOne(p3,p5,p6,M00,M0)        
        .sort

        #call averts(p6,M0)
        .sort

        #call TadpoleM0(x6,p6,M0,0)        
        .sort

* Print+s;
* .end
        
* #include expandnomdeno
        
#endprocedure        



************************************************************
#procedure topm3
*
* topm3
* (simple topology)
*

        #message this is topm3

        if( count(intm3,1)) id p1=p2-p6;

        #call IntOne(p2,p1,p6,m3,M00)                
        .sort

        if( count(intM00,1)) id p4=p3-p6;
        
        #call IntOne(p3,p4,p6,M00,M0)        
        .sort
        #call averts(p6,M0)
        .sort

        #call TadpoleM0(x6,p6,M0,0)        
        .sort

* #include expandnomdeno
        
#endprocedure        


************************************************************
#procedure topm4
*
* topm4
* (simple topology)
*

#message this is topm4
* if( count(intm4,1)) multiply intM00*intM0/intm4;
*         #call one00(p3,p4,p6,e3,e4,e6,int1)

* Numerator p4->p3,p6
        if(count(intm4,1)) id p4=p6-p3;        
        #call IntOne(p3,p4,p6,m4,MxM)
.sort

Print+s;
.sort:AFP-one1;

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


************************************************************
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
        id s1m*s2m*s5m = M^2*miT1*int0/intMMM;

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



************************************************************
#procedure topt1
        #message this is topt1
        
        if( count(intt1,1));
        if( (count(x3,1)=0) && (count(x4,1)=1)
        && (count(x5,1)=1) && (count(x6,1)=1) )
        id x4*x5*x6 =  + (M^2*Gam(-1,1))^3;
        
        Multiply int0/intt1;
        endif;
        
#endprocedure        




************************************************************        
#procedure topn1
        #message this is topn1
        

        
* integrals of type N1
        
        if(count(intn1,1));        
        if( (count(x3,1)=1) && (count(x4,1)=1)
        && (count(x5,1)=1) && (count(x6,1)=1) )
        id x3*x4*x5*x6 =M^4*miBN*int0/intn1;
        
        endif;
        
        .sort 
        
* The following commands should not be used if 'EXACTEXP' is defined.

* id agam=2;
* id bgam=23/3;
* id cgam=35/2+3*z2;
* id dgam=275/192*16+23/2*z2-2*z3;
* id egam=-16*(189/128 - 89/48*z3 - 3/32*z4 - 105/64*z2 - 9/64*z2^2);
* id fgam = -384*(
*          14917/18432 + 1/128*z3*z2 - 175/256*z3 + 649/1536*z4 + 1/320*z5
*          - 275/3072*z2 - 23/1024*z2^2
*          )
*          +16*B4;

* argument;
*   id agam=2;
*   id bgam=23/3;
*   id cgam=35/2+3*z2;
*   id dgam=275/192*16+23/2*z2-2*z3;
*   id egam=-16*(189/128 - 89/48*z3 - 3/32*z4 - 105/64*z2 - 9/64*z2^2);
*   id fgam = -384*(
*          14917/18432 + 1/128*z3*z2 - 175/256*z3 + 649/1536*z4 + 1/320*z5
*          - 275/3072*z2 - 23/1024*z2^2
*          )
*          +16*B4;
* endargument;

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

        #message "                                                     "
        #message "                                                     "
        #message "      _____ _____ _____ _____ ____                   "
        #message "     |     |  _  |_   _|  _  |    \ ___ ___ ___      "
        #message "     | | | |     | | | |     |  |  |___|   | . |     "
        #message "     |_|_|_|__|__| |_| |__|__|____/    |_|_|_  |     "
        #message "                                           |___|     "
        #message "                                                     "
        #message "                                                     "
        #message "     Originally written by M.Steinhauser             "
        #message "     About exact version contact to:                 "
        #message "     pikelner[at]theor.jinr.ru, Andrey Pikelner      "
        #message "                                                     "
        #message MATAD called for `LOOPS'-loop integral

        #call tad`LOOPS'l
        

        if(count(int0,1));
        Multiply 1/int0;  
        
        else;
        exit "Not all integrals reduced";

        endif;        
#endprocedure


#procedure subSimple



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

#endprocedure

#procedure GammaArgToOne
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


#procedure cutep(x)
multiply, 1/ep^('x');
id ep=0;
multiply, ep^('x');
#endprocedure



#procedure subvalues
        
* Expansions from original MATAD        
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
        
* MMM000(1,1,1,1,1,1)        
        id miDM = 2*z3/ep + DM + ep*DMep + ep^3*miDMtrunc;

* MMM000(1,1,1,1,1,0)/M^2        
*** the coefficients of 1/ep^3 and 1/ep^2 are determined
*** from the cancellation of poles 
        id miE3 = -2/3/ep^3-11/3/ep^2 + (-14 + (27*S2)/2 - 2*z2)/ep + E3 + ep*E3ep + ep^3*miE3trunc;

* MM00MM(1,1,0,0,1,1)/M^4                
        id miBN = agam/ep^3+bgam/ep^2+cgam/ep+dgam+egam*ep
        +fgam*ep^2+ggam*ep^3+hgam*ep^4 + ep^5*miBNtrunc;

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
        
* M000MM(1,1,0,0,1,1)
        id miBN1x00 =
        + ep^-3 + ep^-2 * ( 15/4 ) + ep^-1 * ( 65/8 + 12/8*z2 )
        + 81/4*S2 - z3 + 135/16 + 90/16 *z2
        + ep*OepS2
        + ep^2*Oep2S2
        + ep^3*miBN1x00trunc;

* M000MM(1,1,1,1,1,1)
        id miBN1x11 = + (2*z3/ep + D3) + ep*OepD3 + ep^2*miBN1x11trunc;
        

#endprocedure


*--#[ expansion :
*
#procedure expansion(maxeppow)
*         #call subvalues
*
*	Expands the PolyRatFun to sufficient powers in ep.
*
        .sort:expansion-1;
        S DUMMYSYMBOL;        
        PolyRatFun;
* 
* We introduce DUMMYSYMBOL for proper terms ordering in SplitArg
* First term has smallest power of ep        
*         
        id	den(x?) = den(DUMMYSYMBOL*x);
        id	rat(x1?,x2?) = num(x1)*den(DUMMYSYMBOL*x2);
        
        SplitArg,den;
        Multiply replace_(DUMMYSYMBOL,1);
        
**** id	den(?a,x1?) = den(x1,?a);
        repeat id den(x1?,x2?,x3?,?a) = den(x1,x2+x3,?a);
        id	den(x1?,x2?) = den(1,x2/x1)/x1;
        id	den(x1?) = 1/x1;
        
        .sort:expansion-2;
        id	num(x1?) = x1;
*
*       In masters we have minimal power ep^-3
*       And need to expand rat(ep) up to maxpow+3
*         
        if ( count(ep,1) > {`maxeppow' + 3} ) discard;
        repeat;
	        id den(1,x?) = 1-x*den(1,x);
	        if ( count(ep,1) > {`maxeppow'+ 3} ) discard;
        endrepeat;
        .sort:expansion-3;
        Symbol ep(:{`maxeppow'+3});
*         
        #call subvalues

        if ( count(ep,1) > `maxeppow' ) discard;
#endprocedure
*
*--#] expansion : 




*--#[ fillbntab : 
* 
* Example program to fill BNd and BNd0 table
* Can be used if tables for dala^n, n>11
* Inverse of eq.4 used with old routine for
* reduction.        
*         
* #-
* #include matad-ng.hh
* S n3m2;

* * Fill table upto weight 11
* #do dp=0,11
* #do i6=1,2
* #do i5=1,`i6'
* #do i4=1,`i5'
* #do i3=1,`i4'
*         #message dala^`dp'-- `i3' `i4' `i5' `i6'
*         L dalaBN`dp' =dala^`dp'*x3^`i3'*x4^`i4'*x5^`i5'*x6^`i6';
        
*         #do i=1,1
*                 id,once dala^n1?pos_*x3^n3m2?*x4^n4?*x5^n5?*x6^n6? = 
*                 ((4*((n3m2+2)-1)-2*rat(4-2*ep,1))/4/M^2/((n3m2+2)-1)*x3^((n3m2+2)-1)*x4^n4*x5^n5*x6^n6 -
*                 x3^(n3m2+2)*x4^n4*x5^n5*x6^n6)*4*M^2*(n3m2)*((n3m2+2)-1)*dala^(n1-1);
                
*                 if(count(dala,1)) redefine i "0";
*                 .sort
*         #enddo
        
*         Multiply replace_(x3,s1m,x4,s2m,x5,s5m,x6,s6m);
*         #call matad(3)
*         b int0;
*         .sort
        
*         #write<BNtbl> "id,only intbn*dala^`dp'*x3^`i3'*x4^`i4'*x5^`i5'*x6^`i6' = %e\n",dalaBN`dp'
        
*         drop;
*         .sort
* #enddo
* #enddo
* #enddo
* #enddo
* #enddo
* .end
*--#] fillbntab : 
