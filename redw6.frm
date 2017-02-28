#-
#include matad-ng.hh

L exw6 = tad3l([MM0MMM],2,2,2,1,1,1);

* Convert to original MATAD notation
#call FromAuxTopo

* 3 - three-loop reduction
#call matad(3)
.sort

Format mathematica;
#write<redw6.in> "(%E)",exw6
.end
