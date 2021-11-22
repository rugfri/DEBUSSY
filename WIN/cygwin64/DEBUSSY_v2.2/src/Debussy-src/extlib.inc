LIBDIREXT1=/home/tosca01/DEBUSSY_v2.3/extlib/libcerf-1.5/lib/.libs
LIBDIREXT2=/home/tosca01/DEBUSSY_v2.3/extlib/nlopt-2.3/.libs
LIBDIREXT3 =/home/tosca01/DEBUSSY_v2.3/extlib/LAPACK95-lite
INCLDIREXT1 =/home/tosca01/DEBUSSY_v2.3/extlib/LAPACK95-lite/lapack95_modules
####
LIB1 = -L$(LIBDIREXT1) -lcerf
LIB2 = -L$(LIBDIREXT2) -lm -lnlopt
LIB3 = -L$(LIBDIREXT3) -llapack95 -llapack -lblas
