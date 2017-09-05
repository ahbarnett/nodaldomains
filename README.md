# nodaldomains

extract nodal domains, their genus and percolation properties, from 2D or 3D uniformly gridded functions. C++ code with MATLAB/octave wrappers.

Main author: Alex Barnett (c) 2017

Code contributions/acknowledgments: Kyle Konrad (2011-2012), Matthew Jin (2014-2015), Ziff-Newman, Tarjan, Sedgewick percolation algorithms (2001).

## Dependencies

- C++ compiler and GNU make
- If you want interfaces to them: [MATLAB](http://mathworks.com) and/or octave (including its development libraries)

## Installation

From the main directory:

`make test` compiles the library `domainlib.o` then tests it via C++ drivers.

`make octave` builds and tests the octave interfaces.

`make matlab` builds the MATLAB interface; use `test_nodal3dziff.m` and `test_perc3d.m` to test.

Run `make` to see other options.

## Codes

The following routines are available in the library, and MATLAB/octave functions:

- `nodal3dziff` : returns labeled nodal domains, and list of their sizes, for a real function given on 3D cubical grid, using quick union/find.

- `perc3d` : estimates percolation threshold h0 for real function given on 3D cubical grid, and nodal domain information at threshold, using quick union/find in style of Newman-Ziff.

- `genus3djin` : extracts genus and neighbor relations from labeled nodal domains on a 3D cubical grid.
