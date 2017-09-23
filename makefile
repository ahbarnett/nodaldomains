# makefile for nodal statistics utilities
# Barnett 8/30/17, added 2d 9/22/17

CXX=g++
CXXFLAGS = -fPIC -Ofast -funroll-loops -march=native

# location of MATLAB's mex compiler...
MEX=mex
# Mac users should use something like this:
#MEX = /Applications/MATLAB_R2017a.app/bin/mex
# location of your MWrap executable...
MWRAP=mwrap

.PHONY: test clean mexclean all octave matlab

default: usage

usage:
	@echo "Makefile for nodal domains library. Specify what to make:"
	@echo "make lib - compile domainlib.o"
	@echo "make test - compile lib and run test drivers in C++"
	@echo "make matlab - compile MEX interface"
	@echo "make octave - compile and test octave interface"
	@echo "make all - does the above 3 things"
	@echo "make clean - remove objects and executables"
	@echo "make mex - recreate MEX interface from .mw (needs mwrap)"
	@echo "make mexclean - remove generated MEX interface (experts only)"

%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

lib: domainlib.o

test: driver2d driver3d
	./driver2d 1
	./driver3d 1

all: test matlab octave

driver2d: driver2d.cpp domainlib.o
	$(CXX) $(CXXFLAGS) driver2d.cpp -o driver2d domainlib.o -lm

driver3d: driver3d.cpp domainlib.o
	$(CXX) $(CXXFLAGS) driver3d.cpp -o driver3d domainlib.o -lm

matlab: gateway.cpp domainlib.h domainlib.o
	$(MEX) gateway.cpp domainlib.o -largeArrayDims -lm

octave: gateway.cpp domainlib.h domainlib.o
	mkoctfile --mex gateway.cpp domainlib.o -lm
	octave test_nodal2dziff.m
	octave test_nodal3dziff.m
	octave test_perc2d.m
	octave test_perc3d.m

mex: domainlib.mw
	$(MWRAP) -list -mex gateway -mb domainlib.mw
	$(MWRAP) -mex gateway -c gateway.cpp domainlib.mw

clean:
	rm -f *.o driver2d driver3d
# only do this if you have mwrap to rebuild the interfaces...
mexclean: clean
	rm -f gateway.cpp nodal2dziff.m nodal3dziff.m perc2d.m perc3d.m gateway.mex*
