MATAD-ng
====

**MATAD** is a form code for three-loop MAssive TADpole integral reduction originally written by M.Steinhauser and published in [Comput.Phys.Commun. 134 (2001) 335-364](http://inspirehep.net/record/532857).

**MATAD-ng** is based on the same code by M.Steinhauser but use FORM4 PolyRatFun for intermediate expressions and allow to produce reduction without expansion in ep in terms of master integrals and Euler gamma functions.

![MATAD-ng topologies](https://raw.githubusercontent.com/wiki/apik/matad-ng/images/topmtd.png)

For origianl MATAD usage and examples see M.Steinhauser lectures at CAPP 2009 https://www.ttp.kit.edu/~ms/capp09.pdf

To use reduction at arbitrary **d** only following 9 master integrals needed:
![MATAD-ng master integrals](https://raw.githubusercontent.com/wiki/apik/matad-ng/images/masterints.png)

With expansion of Euler gamma functions near d for remaining trivial integrals

> `Gam(n,x)=Gamma(n+(2-d/2)*x)*Exp(ep*x*EulerGamma), iGam(n,x)=Exp(-ep*x*EulerGamma)/Gamma(n+(2-d/2)*x)`

Each integral is divided by `Exp(-ep*EulerGamma)` and `Zeta[2]` is present in final results.
