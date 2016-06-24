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


* Input in terms of origianl MATAD variables s1m,s2m,s3m,p1,p2,p3
L ex2loopA = s1m^2*s2m*s3m*p2.p3;

* And in using auxiliary topology 
L ex2loopB = tad2l([MM0],2,2,2);

* Convert to original MATAD notation
#call FromAuxTopo

* 2 - two-loop reduction
#call matad(2)
* expansion upto ep^2
#call exp4d(2)
hide;
.sort



* Input in terms of origianl MATAD variables s1m...s6m,p1...p6
L ex3loopA = s1m^2*s2m*s3m*p2.p3*s4m*s5m*p6.p6;

* And in using auxiliary topology 
L ex3loopB = tad3l([MM00MM],2,2,2,1,1,1);

* Convert to original MATAD notation
#call FromAuxTopo

* 3 - three-loop reduction
#call matad(3)
* expansion upto ep^1
#call exp4d(1)
hide;
.sort

* Now print results:

unhide;

b ep,M;
Print+s;
.end
