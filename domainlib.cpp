#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "domainlib.h"

using namespace std;

#define SIGN(x) (x>=0.0)

int findroot(int* ptr, int i);
void qunion(int* ptr, int p, int q);

int findroot(int* ptr, int i)
// Returns the label of the root of the tree that site i is in.
// Non-recursive, Tarjan path halving version, by Ziff 2001, without globals
{
  int r,s;
  r = s = i;
  while (ptr[r]>=0) {  // negative ptr indicates it's a root
    ptr[s] = ptr[r];
    s = r;
    r = ptr[r];
  }
  return r;
}

void qunion(int* ptr, int p, int q)
// Joins two clusters (trees), namely p and q, using the pointer array ptr.
// Adapted from quick-union in Sedgewick notes on union/find. Barnett.
{
  int i = findroot(ptr, p);
  int j = findroot(ptr, q);
  if (i==j) return;         // nothing to do
  if (-ptr[i] < -ptr[j]) {
    ptr[j] += ptr[i];       // add the sizes (both negative)
    ptr[i] = j;             // make i child of q's root.
  } else {                  // do it other way around...
    ptr[i] += ptr[j];
    ptr[j] = i;
  }
}

int nodal3dziff(int N, double *u, int *d, int *siz, int *nd, int verb)
/* Simple Ziff quick union/find alg for labeling +/- nodal domains on 3D grid.
   Single thread only.

   Inputs: u, a real NxNxN contiguous array of doubles.
           N, the cube side length (positive integer).
   Output: d, an int NxNxN contiguous array containing domain numbers
              where domains are regions in u with same sign, connected along
              the usual six grid point cube faces. Numbers are in 1,..,nd.
           siz, a nd-length list of sizes of the domains.
           nd, the number of domains = max(d).

   Barnett 8/30/17
*/
{
  int n=N*N*N;    // # grid pts
  int *ptr;
  ptr = (int*)malloc((sizeof(int))*n);
  if (verb) printf("nodal3dziff: N=%d, n=%d\n",N,n);
  
  //for (int i=0;i<n;++i) d[i] = (u[i] >= 0) ? 1 : -1; // sign func, as warmup
  
  for (int i=0;i<n;++i) ptr[i] = -1;       // each site set to -(cluster size)

  int i=0;    // grid pts (sites)
  for (int z=0; z<N; ++z)
    for (int y=0; y<N; ++y)
      for (int x=0; x<N; ++x) {
	int s = SIGN(u[i]);    // site sign
	// use connectedness given by signs, without periodic wrapping:
	if (x>0) {
	  int j = i-1;         // x neighbor
	  if (SIGN(u[j])==s) qunion(ptr,i,j);
	}
	if (x<N-1) {
	  int j = i+1;         // x neighbor
	  if (SIGN(u[j])==s) qunion(ptr,i,j);
	}
	if (y>0) {
	  int j = i-N;         // y neighbor
	  if (SIGN(u[j])==s) qunion(ptr,i,j);
	}
	if (y<N-1) {
	  int j = i+N;         // y neighbor
	  if (SIGN(u[j])==s) qunion(ptr,i,j);
	}
	if (z>0) {
	  int j = i-N*N;         // z neighbor
	  if (SIGN(u[j])==s) qunion(ptr,i,j);
	}
	if (z<N-1) {
	  int j = i+N*N;         // z neighbor
	  if (SIGN(u[j])==s) qunion(ptr,i,j);
	}
       	++i;
      }
  if (verb) printf("unions done\n");
  
  int *dn;
  dn = (int*)malloc(sizeof(int)*n);     // domain numbers, need for all sites!
  int dc=0;                             // domain counter
  for (int j=0; j<n; ++j) {
    if (ptr[j]<0) {                     // only write to the root sites
      siz[dc] = -ptr[j];
      dn[j] = dc+1;                     // 1-indexed
      if (verb>2 && dc<50) printf("dn[%d] = %d, siz = %d\n",j,dc+1,siz[dc]); // debug
      dc++;
    }
  }
  *nd = dc;                            // write to output arg
  
  for (int i=0; i<n; ++i)              // find all sites' domain labels...
    d[i] = dn[findroot(ptr,i)];        //  ...mapping to their 1,..,nd numbers.

  if (verb) printf("domains labels done (nd=%d)\n",*nd);
  
  free(ptr);
  free(dn);
  return 0;
}


