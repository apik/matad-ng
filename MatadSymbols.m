(*  

   [1] MATAD: A Program package for the computation of MAssive TADpoles 
       Comput.Phys.Commun. 134 (2001) 335-364, hep-ph/0009029 
   [2] High-precision epsilon expansions of single-mass-scale four-loop vacuum bubbles
       JHEP 0506 (2005) 051, hep-ph/0503209
   [3] Relations for Nielsen polylogarithms
       J.Approx.Theor. 193 (2015) 74-88

*)

(* Functions from [3] *)
NielsenQ[L_]:=Length[L]>=1&&Complement[Rest[L],{1}]==={};
N[Li[L_List,x_],prec_:$MachinePrecision]:=N[PolyLog[L[[1]]-1,Length[L],x],prec]/;NielsenQ[L];
N[Cl[L_List,x_],prec_:$MachinePrecision]:=Module[{p=N[PolyLog[L[[1]]-1,Length[L],Exp[I x]],prec]},If[EvenQ[Total[L]],Im[p],Re[p]]]/;NielsenQ[L];
N[Gl[L_List,x_],prec_:$MachinePrecision]:=Module[{p=N[PolyLog[L[[1]]-1,Length[L],Exp[I x]],prec]},If[EvenQ[Total[L]],Re[p],Im[p]]]/;NielsenQ[L];

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

        DN -> 6*Zeta[3] - 4*Zeta[2]*Log[2]^2 + 2/3Log[2]^4 - 21/2*Zeta[4] - 16*Li[{2},1/2],
        
        B4 -> -4*Zeta[2]Log[2]^2 + 2/3*Log[2]^4 - 13/2*Zeta[4] + 16*Li[{4},1/2],

        (* eq:6.23 from [2] *)
        E3 -> -139/3 - (7*Pi^2)/4 - (4*Pi^3)/(9*Sqrt[3]) + 15*Sqrt[3]*Cl[{2}, (2*Pi)/3] -
        6*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3] - 3*Sqrt[3]*Cl[{2}, (2*Pi)/3]*Log[3] + Zeta[3]/3,
        
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
        6*Sqrt[3]*Gl[{2, 1}, (2*Pi)/3] - 3*Sqrt[3]*Cl[{2}, (2*Pi)/3]*Log[3] + Zeta[3]
    }


miSubs=
    {
        miT1 ->
        -21/2 - 3/(2*ep^2) - 9/(2*ep) + (27*S2)/2 - (3*z2)/2 
        + ep * T1ep + ep^2 * T1ep2 + ep^3*Trunc[miT1],
        
        miD6 -> 2*z3/ep + D6 + ep*D6ep + ep^2*Trunc[miD6],
        
        miD5 -> 2*z3/ep + D5 + ep*D5ep + ep^2*Trunc[miD5],
        
        miD4 -> 2*z3/ep + D4 + ep*D4ep + ep^2*Trunc[miD4],
        
        miDN -> 2*z3/ep + DN + ep*DNep + ep^2*Trunc[miDN],
        
        miDM -> 2*z3/ep + DM + ep*DMep + ep^2*Trunc[miDM],

        miE3 -> -2/3/ep^3-11/3/ep^2 + (-14 + (27*S2)/2 - 2*z2)/ep + E3 + ep*E3ep + ep^2*Trunc[miE3],

        miBN -> 275/12 + 2/ep^3 + 23/(3*ep^2) + ep^3*ggam + ep^4*hgam + ep^5*Trunc[miBN] + (23*z2)/2 +
        (35/2 + 3*z2)/ep - 2*z3 - 16*ep*(189/128 - (105*z2)/64 - (9*z2^2)/64 - (89*z3)/48 - (3*z4)/32) +
        ep^2*(16*B4 - 384*(14917/18432 - (275*z2)/3072 - (23*z2^2)/1024 - (175*z3)/256 + (z2*z3)/128 + (649*z4)/1536 + z5/320)),
                
        miBN1x00 ->
        + ep^-3 + ep^-2 * ( 15/4 ) + ep^-1 * ( 65/8 + 12/8*z2 )
        + 81/4*S2 - z3 + 135/16 + 90/16 *z2
        + ep*OepS2
        + ep^2*Oep2S2
        + ep^3*Trunc[miBN1x00],

        miBN1x11 -> + (2*z3/ep + D3) + ep*OepD3 + ep^2*Trunc[miBN1x11]
    }