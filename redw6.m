(* Load reasults of master inegrals epsilon expansion *)
Get["mtdw6.m"];
(* HPL's with 20k digits accuracy *)
nhpl=Join[Get["nhplRe.m"],Get["nhplIm.m"]];
(* Result of reduction *)
expr=Get["redw6.in"];

(* Expansion *)
exprExp=Collect[FunctionExpand[Normal[Series[expr//.{d->4-2ep,rat[a_,b_]->(a/b)}/.miw6,{ep,0,1}]]],{M,ep},Expand];

Print["Befor numerical evaluation: ",InputForm[exprExp]];

exprNum = N[exprExp/.subuniweight/.nhpl,100];

Print["Numerical value: ",InputForm[exprNum]];
