(*  

   [1] MATAD: A Program package for the computation of MAssive TADpoles 
       Comput.Phys.Commun. 134 (2001) 335-364, hep-ph/0009029 
   [2] High-precision epsilon expansions of single-mass-scale four-loop vacuum bubbles
       JHEP 0506 (2005) 051, hep-ph/0503209
   [3] Relations for Nielsen polylogarithms
       J.Approx.Theor. 193 (2015) 74-88
   [4] High precision calculation of multiloop Feynman integrals by difference equations
       Int.J.Mod.Phys. A15 (2000) 5087-5159, hep-ph/0102033
   [5] New results for the epsilon expansion of certain one, two and three loop Feynman diagrams
       Nucl.Phys. B605 (2001) 266-318, hep-th/0012189 

*)

(* Functions from [3] *)
NielsenQ[L_]:=Length[L]>=1&&Complement[Rest[L],{1}]==={};
N[Li[L_List,x_],prec_:$MachinePrecision]:=N[PolyLog[L[[1]]-1,Length[L],x],prec]/;NielsenQ[L];
N[Cl[L_List,x_],prec_:$MachinePrecision]:=Module[{p=N[PolyLog[L[[1]]-1,Length[L],Exp[I x]],prec]},If[EvenQ[Total[L]],Im[p],Re[p]]]/;NielsenQ[L];
N[Gl[L_List,x_],prec_:$MachinePrecision]:=Module[{p=N[PolyLog[L[[1]]-1,Length[L],Exp[I x]],prec]},If[EvenQ[Total[L]],Re[p],Im[p]]]/;NielsenQ[L];

miSubs=
    {
        miT1 ->
        -21/2 - 3/(2*ep^2) - 9/(2*ep) + (27*S2)/2 - (3*Zeta[2])/2 
        + ep * T1ep + ep^2 * T1ep2 + ep^3 * T1ep3 + ep^4*Trunc[T1ep4],
        
        miD6 -> 2*Zeta[3]/ep + D6 + ep*D6ep + ep^2*Trunc[D6ep2],
        
        miD5 -> 2*Zeta[3]/ep + D5 + ep*D5ep + ep^2*Trunc[D5ep2],
        
        miD4 -> 2*Zeta[3]/ep + D4 + ep*D4ep + ep^2*Trunc[D4ep2],
        
        miDN -> 2*Zeta[3]/ep + DN + ep*DNep + ep^2*Trunc[DNep2],
        
        miDM -> 2*Zeta[3]/ep + DM + ep*DMep + ep^2*Trunc[DMep2],

        miE3 -> -2/3/ep^3-11/3/ep^2 + (-14 + (27*S2)/2 - 2*Zeta[2])/ep + E3 + ep*E3ep + ep^2*Trunc[E3ep2],

        miBN ->  2/ep^3 + 23/(3*ep^2) + (35/2 + 3*Zeta[2])/ep + 275/12 + (23*Zeta[2])/2  - 2*Zeta[3] -
        16*ep*(189/128 - (105*Zeta[2])/64 - (9*Zeta[2]^2)/64 - (89*Zeta[3])/48 - (3*Zeta[4])/32) +
        ep^2*(16*B4 - 384*(14917/18432 - (275*Zeta[2])/3072 - (23*Zeta[2]^2)/1024 - (175*Zeta[3])/256 +
                           (Zeta[2]*Zeta[3])/128 + (649*Zeta[4])/1536 + Zeta[5]/320)) + ep^3*BNep3 + ep^4*BNep4 + ep^5*Trunc[BNep5],
                
        miBN1x00 ->
        + ep^-3 + ep^-2 * ( 15/4 ) + ep^-1 * ( 65/8 + 12/8*Zeta[2] )
        + 81/4*S2 - Zeta[3] + 135/16 + 90/16 *Zeta[2]
        + ep*OepS2
        + ep^2*Oep2S2
        + ep^3*Trunc[BN1x00ep3],

        miBN1x11 -> + (2*Zeta[3]/ep + D3) + ep*OepD3 + ep^2*Trunc[BN1x11ep2],

        Gam[x_, y_]  -> Exp[ep*y*EulerGamma]*Gamma[x + ep*y],
        iGam[x_, y_] -> Exp[-ep*y*EulerGamma]/Gamma[x + ep*y] 
    };

symbSubs =
    {
        D6 -> 6*Zeta[3] - 17*Zeta[4] - 4*Zeta[2]*Log[2]^2 + 2/3*Log[2]^4 + 16*PolyLog[4,1/2] -
        4*Cl[{2},Pi/3]^2,

        (* eq:6.25 from [2] *)
        D5 -> (-4*Pi^4)/405 - 6*Cl[{2}, (2*Pi)/3]^2 + 12*Gl[{3, 1}, (2*Pi)/3] + 6*Zeta[3] -
        (26*Log[3]*Zeta[3])/3,
        
        D4 -> 6*Zeta[3] - 77/12*Zeta[4] - 6*Cl[{2},Pi/3]^2,

        D3 -> 6*Zeta[3] - 15/4*Zeta[4] - 6*Cl[{2},Pi/3]^2,

        DM -> 6*Zeta[3] - 11/2*Zeta[4] - 4*Cl[{2},Pi/3]^2,

        DN -> 6*Zeta[3] - 4*Zeta[2]*Log[2]^2 + 2/3*Log[2]^4 - 21/2*Zeta[4] + 16*Li[{4},1/2],
        
        B4 -> -4*Zeta[2]Log[2]^2 + 2/3*Log[2]^4 - 13/2*Zeta[4] + 16*Li[{4},1/2],

        (* eq:6.13 from [2] *)
        BNep3 -> -48005/32 - (189*Pi^2)/32 - (4219*Pi^4)/192 + (631*Pi^6)/60480 + (272*Pi^4*Log[2])/15 -
        80*Pi^2*Log[2]^2 + (64*Pi^2*Log[2]^3)/3 + 80*Log[2]^4 - (64*Log[2]^5)/5 +
        1920*PolyLog[4, 1/2] + 1536*PolyLog[5, 1/2] + (15965*Zeta[3])/12 +
        (89*Pi^2*Zeta[3])/12 + Zeta[3]^2 - (6223*Zeta[5])/5,
        
        (* eq:6.23 from [2] *)
        E3 -> -139/3 - (7*Pi^2)/4 - (4*Pi^3)/(9*Sqrt[3]) + 15*Sqrt[3]*Cl[{2}, (2*Pi)/3] -
        6*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3] - 3*Sqrt[3]*Cl[{2}, (2*Pi)/3]*Log[3] + Zeta[3]/3,

        (* eq:4.25 from [5], O(ep) part of eq:6.23 from [2] is incorrect *)
        E3ep -> -430/3 - (19*Pi^2)/3 - (20*Pi^3)/(9*Sqrt[3]) - (37*Pi^4)/720 + 
        34*Sqrt[3]*Cl[{2}, Pi/3] + (13*Pi^2*Cl[{2}, Pi/3])/(6*Sqrt[3]) + 
        (Pi^2*Cl[{2}, (2*Pi)/3])/(2*Sqrt[3]) + (40*Cl[{4}, Pi/3])/Sqrt[3] - 
        3*Sqrt[3]*Cl[{4}, (2*Pi)/3] + 12*Sqrt[3]*Cl[{2, 1, 1}, (2*Pi)/3] - 
        30*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3] + (4*Pi^3*Log[3])/(9*Sqrt[3]) - 
        10*Sqrt[3]*Cl[{2}, Pi/3]*Log[3] + 6*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3]*Log[3] + 
        Sqrt[3]*Cl[{2}, Pi/3]*Log[3]^2 - (2*Zeta[3])/3 - (14*Pi*Zeta[3])/(3*Sqrt[3]) + 
        3*Sqrt[3]*Pi*Zeta[3],
        
        S2 -> 4*Cl[{2},Pi/3]/9/Sqrt[3],

        (* eq:6.12 from [2] *)
        OepS2 -> -763/32 + (Pi^2*(975 - 160*Sqrt[3]*Pi + 19*Pi^2))/480 -
        27*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3] - (27*Sqrt[3]*Cl[{2}, (2*Pi)/3]*(-5 + Log[9]))/4 -
        (15*Zeta[3])/4,

        (* eq:6.12 from [2] *)
        Oep2S2 -> -36*Sqrt[3]*Cl[{4}, Pi/3] - (81*Sqrt[3]*Cl[{4}, (2*Pi)/3])/2 +
        162*Sqrt[3]*Cl[{2, 1, 1}, (2*Pi)/3] +
        (9*Sqrt[3]*Cl[{2}, (2*Pi)/3]*(145 + 3*Pi^2 + 18*(-5 + Log[3])*Log[3]))/8 +
        (81*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3]*(-5 + Log[9]))/2 +
        (-127050 + 1350*Pi^2 + 95*Pi^4 + 320*Sqrt[3]*Pi^3*(-5 + Log[9]) +
         80*(-65 + 156*Sqrt[3]*Pi - 2*Pi^2)*Zeta[3] - 384*Zeta[5])/640,

        (* eq:6.9 from [2] *)
        T1ep -> -45/2 - (3*Pi^2)/4 - (2*Pi^3)/(9*Sqrt[3]) + 9*Sqrt[3]*Cl[{2}, (2*Pi)/3] -
        6*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3] - 3*Sqrt[3]*Cl[{2}, (2*Pi)/3]*Log[3] + Zeta[3],

        T1ep2 -> -93/2 - (7*Pi^2)/4 - (7*Pi^4)/240 - 3*Sqrt[3]*Cl[{4}, (2*Pi)/3] +
        12*Sqrt[3]*Cl[{2, 1, 1}, (2*Pi)/3] + (2*Pi^3*(-3 + Log[3]))/(9*Sqrt[3]) +
        6*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3]*(-3 + Log[3]) +
        (Cl[{2}, (2*Pi)/3]*(4*Pi^2 + 9*(14 + (-6 + Log[3])*Log[3])))/(2*Sqrt[3]) +
        (3 + (13*Pi)/(3*Sqrt[3]))*Zeta[3],
        
        T1ep3 -> -189/2 - (7*Pi^4)/80 - (187*Pi^5)/(1620*Sqrt[3]) + 6*Sqrt[3]*Gl[{4, 1}, (2*Pi)/3] -
        24*Sqrt[3]*Gl[{2, 1, 1, 1}, (2*Pi)/3] + 3*Sqrt[3]*Cl[{4}, (2*Pi)/3]*(-3 + Log[3]) -
        12*Sqrt[3]*Cl[{2, 1, 1}, (2*Pi)/3]*(-3 + Log[3]) - (Pi^3*(14 + (-6 + Log[3])*Log[3]))/
        (9*Sqrt[3]) - 3*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3]*(14 + (-6 + Log[3])*Log[3]) -
        (Pi^2*(45 + 16*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3] - 2*Zeta[3]))/12 + 7*Zeta[3] +
        (18*Pi*Gl[{3, 1}, (2*Pi)/3] - 13*Pi*(-3 + Log[3])*Zeta[3])/(3*Sqrt[3]) +
        (Cl[{2}, (2*Pi)/3]*(-4*Pi^2*(-3 + Log[3]) - 3*(-90 + Log[3]*(42 + (-9 + Log[3])*Log[3]) + 4*Zeta[3])))/(2*Sqrt[3]) +
         (3*Zeta[5])/5
    };


(* Numerical values for higher orders in ep *)

N[D6ep,prec_:$MachinePrecision] = 41.87670208303`13;
N[D5ep,prec_:$MachinePrecision] = 74.9083038310`13;
(* eq:4.4 from [2] *)
N[D4ep,prec_:$MachinePrecision]  = 31.79387586520335091502703130598293231890424270640`50;
N[D4ep2,prec_:$MachinePrecision] = -95.53186858506048154153099646099149596350732620450`50;
(* eq:4.6 from [2] *)
N[DNep,prec_:$MachinePrecision]  = 30.6810352753459`15;
N[DNep2,prec_:$MachinePrecision] = -13.303460640856`15;
(* eq:4.5 from [2] *)
N[DMep,prec_:$MachinePrecision]  = 29.00667443783775908331902631581722417518004677385`50;
N[DMep2,prec_:$MachinePrecision] = -62.36139634229625160648107039345983094078357810702`50;
(* eq:4.3 from [2] *)
N[E3ep2,prec_:$MachinePrecision] = -574.65405761296725286615112197885868312967183890`50;
(* eq:198 from [4] *)
N[BNep4,prec_:$MachinePrecision] = -4108.81596022`13;
