#-
#include matad-ng.hh

L exShift1 = tad3l([MMMMMM],2,2,2,1,1,1);

* Convert to original MATAD notation
#call FromAuxTopo

* 3 - three-loop reduction
#call matad(3)

* shift dimension of integral to (4-2) - 2*ep and rewrite
* it in terms of four dimensional integrals
#call shift4plus(-2)

* expansion upto ep^1
#call exp4d(1)

b ep,M;
Print+s;
.store



L exShift2 = tad3l([MMMMMM],2,2,2,1,1,1);

* Convert to original MATAD notation
#call FromAuxTopo

* 3 - three-loop reduction
#call matad(3)
.sort

* Save reduction valid for arbitrary D in ex4d
L ex4d = exShift2;
.sort

hide ex4d;
.sort

* shift dimension of integral to (4-2) - 2*ep and rewrite
* it in terms of four dimensional integrals
#call shift4plus(-2)

* Shift dimension of all integrals back
#call shift4plus(2)

.sort
unhide;
* expansion up to ep^1
#call exp4d(1)
.sort

* Check that difference is zero
L diffShift = ex4d - exShift2;

b ep,M;
Print+s;
.end
