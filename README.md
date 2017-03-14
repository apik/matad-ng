MATAD-ng
====

**MATAD** is a form code for three-loop MAssive TADpole integral reduction originally written by M.Steinhauser and published in [Comput.Phys.Commun. 134 (2001) 335-364](http://inspirehep.net/record/532857).

**MATAD-ng** is based on the same code by M.Steinhauser but use FORM4 PolyRatFun for intermediate expressions and allow to produce reduction without expansion in ep in terms of master integrals and gamma functions.

![MATAD-ng topologies](https://raw.githubusercontent.com/wiki/apik/matad-ng/images/topmtd.png)

For origianl MATAD usage and examples see M.Steinhauser lectures at CAPP 2009 https://www.ttp.kit.edu/~ms/capp09.pdf

## Master integrals

To use reduction at arbitrary **d** only following 9 master integrals needed:

![MATAD-ng master integrals](https://raw.githubusercontent.com/wiki/apik/matad-ng/images/masterints.png)

With expansion of Euler gamma functions near d for remaining trivial integrals

> `Gam(n,x)=Gamma(n+(2-d/2)*x)*Exp(ep*x*EulerGamma), iGam(n,x)=Exp(-ep*x*EulerGamma)/Gamma(n+(2-d/2)*x)`

Each integral is divided by `Exp(-ep*EulerGamma)` and `Zeta[2]` is present in final results.

## Extension to weight 6

To use new results for integrals expansion up to weight 6 disable expansion in form and instead use attached mathematica files with master integrals substitution rules `mtdw6.m` and tables with HPL's of up to weight 6 of argument `Exp[I Pi/3]`.
Real parts stored in `nhplRe.m` and imaginary in `nhplIm.m` with 20000 numerical digits precision.

## Usage

### Reduction

Calculate one-loop integral using MATAD notation for massive propagators **s1m** and for massles propagator **1/p1.p1**
```
#-
#include matad-ng.hh
* We use s1m for massive denominator and 1/p1.p1 for massless
L ex1loop = s1m^2*p1.p1;
* 1 - one-loop reduction
#call matad(1)
* expansion upto ep^3
#call exp4d(3)
hide;
.sort
```

Calculate three-loop integral using auxiliary topology definition, mass distribution determined by symbol for example `[MM00MM]`

```
L ex3loopB = tad3l([MM00MM],2,2,2,1,1,1);
* Convert to original MATAD notation
#call FromAuxTopo
* 3 - three-loop reduction
#call matad(3)
```

Expand result in **ep** up to O(ep) near d=4-2e*p
```
#call exp4d(1)
```

### Dimension shifts

Here we express two-loop master integral **miT1** in `d=(4+2)-2*ep` dimensions in terms of integrals in `d=4-2*ep` dimensions
```
L f1 = miT1;
#call shift4plus(2)
```
And result is
```
f1 =
   + Gam(1,1)^2*rat(-48,d^6 - 13*d^5 + 64*d^4 - 148*d^3 + 160*d^2 - 64*d)
   + miT1*rat(3,d^2 - d);
```
where **miT1** is two-loop master integral in four dimensions and `d=4-2*ep`

## Tests

To check correctness of results with package supplied test cases with results obtained with the help of original MATAD package. To run tests use command
`$ form tests.frm` and file tests.err will be produced if results of calculation do not match.
