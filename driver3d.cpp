// driver/tester for nodal3dziff, using iid bond percolation on cube of Z^3.
// See makefile for compilation.
// usage: ./driver [verb]
// where verb = 0,1,2,3 sets verbosity of text output.
// Barnett 8/30/17

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdlib.h>
#include "domainlib.h"

using namespace std;

// Random numbers: crappy unif random number generator in [0,1):
#define rand01() ((double)rand()/RAND_MAX)
// unif[-1,1]:
#define randm11() (2*rand01() - (double)1.0)

int main(int argc, char* argv[]) {
  int N=200;
  int verb = 1;
  if (argc>1)
    sscanf(argv[1],"%d",&verb);
  int n=N*N*N;
  vector<double> u(n);
  vector<int> d(n);
  vector<int> siz(n);          // alloc for at most n domains
  int nd;

  srand(time(NULL));           // new seed each run
  for (int i=0;i<n;++i)
    u[i] = rand01() - 0.5;     // signs will give site percolating Z^3 lattice

  // test nodal3dziff :
  int sign = 1;
  int ier = nodal3dziff(N, &u[0], &d[0], &siz[0], &nd, sign, verb);
  if (verb>1) {
    printf("first few entries of u and d...\n");
    for (int i=0;i<50;++i)
      cout << u[i] << "   \t" << d[i] << endl;
    printf("first few i, siz[i] pairs... (i shown 1-indexed to match d)\n");
    if (nd>20) nd=20;         // how many to show; should be 1 or 2 big ones
    for (int i=0;i<nd;++i)
      cout << i+1 << "    \t" << siz[i] << endl;
  }

  // test perc3d :
  int nh = 100000;                 // other than malloc, no CPU cost to large
  double h0, hran[2] = {0.0,0.45};  // h interval to sweep in nh steps
  ier = perc3d(N, &u[0], &d[0], &siz[0], &nd, nh, hran, &h0, verb);
  if (verb) {
    printf("h0 = %g\n",h0);
    if (h0==hran[1])
      printf("h0 at top of range: perc threshold probably >h0\n");
  }

  return 0;
}
