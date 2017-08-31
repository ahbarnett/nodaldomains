# nodaldomains

extract nodal domains, their genus and percolation properties, from 2D or 3D uniformly gridded functions. C++ code with MATLAB/octave wrappers.

Main author: Alex Barnett

Code contributions/acknowledgments: Kyle Konrad (2011-2012), Matthew Jin (2014-2015), Ziff-Newman percolation codes (2001).

## Dependencies

- C++ compiler and GNU make
- If you want interfaces to them: [MATLAB](http://mathworks.com) and/or octave

## Installation

`make test` compiles the library `domainlib.o` and tests it via C++ drivers.

`make matlab; make octave` builds the interfaces.

Run `make` to see other options.

## Codes

The following routines are available in the library, and MATLAB/octave functions:

- `nodal3dziff` : returns labeled nodal domains, and list of their sizes, for a real function given on 3D cubical grid, using quick union/find.

