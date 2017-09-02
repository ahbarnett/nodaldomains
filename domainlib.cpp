#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
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


int labeldoms(int n, int* ptr, double *u, int* d, int* siz, int* nd, int sign, double h, int verb)
// Utility to label domains with numbers 1,..,nd. Doesn't care about dimension.
// If sign!=0, omits domains whose root has u<h (if sign<0) or u>h (if sign>0).
// Outputs: d, grid of domain numbers; siz, domain sizes; nd # of domains.
{
  int *dn;
  dn = (int*)malloc(sizeof(int)*n);     // domain numbers, need for all sites!
  int dc=0;                             // domain counter
  for (int j=0; j<n; ++j) {
    if (ptr[j]<0)                     // only write to the root sites
      if (sign>0 && u[j]<h || sign<0 && u[j]>h)
	dn[j] = 0;                     // not a domain; don't increment dc
      else {
	siz[dc] = -ptr[j];
	dn[j] = dc+1;                     // 1-indexed
	// for debug only...
	if (verb>2 && dc<50) printf("dn[%d] = %d, siz = %d\n",j,dc+1,siz[dc]);
	dc++;
      }
  }
  *nd = dc;                            // write to output arg
  
  for (int i=0; i<n; ++i)              // find all sites' domain labels...
    d[i] = dn[findroot(ptr,i)];        //  ...mapping to their 1,..,nd numbers.

  if (verb) printf("domains labels done (nd=%d)\n",*nd);
  free(dn);
}

int nodal3dziff(int N, double *u, int *d, int *siz, int *nd, int sign, int verb)
/* Simple Ziff quick union/find alg for labeling +/- nodal domains on 3D grid.
   Single thread only.

   Inputs: u, a real NxNxN contiguous array of doubles.
           N, the cube side length (positive integer).
	   sign: if 0 count domains of both signs, if >0 just +ve, <0 just -ve.
	   verb: verbosity (0,1,2,3).
   Outputs:
           d, an int NxNxN contiguous array containing domain numbers
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
  for (int z=0; z<N; ++z)  // cluster all nei w/ same sign (ignore sign arg)
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
  
  int ier = labeldoms(n, ptr, u, d, siz, nd, sign, 0.0, verb);
  free(ptr);
  return ier;
}

int binsort(int n, double* u, int nb, double blo, double iwid, vector<int>& inds, vector<int>& offs, int verb)
/* bin-sort a length-n array u using nb bins of 1/width iwid starting at blo,
   and a single bin [blo+nb/iwid, +infty).
   Outputs:
     inds = index list (length = total counts in all bins <= n)
     offs = bin offset list (length nb+2), with entry offs[nb] = total counts in
            the finite bins, and offs[nb+1] this total including up to +infty.
*/
{
  double bhi = blo + nb/iwid;              // top of last finite bin
  vector<int> counts(nb+1,0);              // first find counts in each bin...
  for (int i=0;i<n;++i)
    if (u[i]>blo) {
      int b = (u[i]>bhi) ? nb : (int)floor(iwid*(u[i]-blo)); // bin # for u[i]
      counts[b]++;
    }
  
  offs[0] = 0;
  for (int b=0;b<nb+1;++b)                 // offsets = cumulative sum of counts
    offs[b+1] = offs[b]+counts[b];         // note [nb] entry is total counts
  vector<int> inv(n,-1);                   // inverse reordering (-1 = empty)
  vector<int> tempoffs = offs;             // array copy
  for (int i=0;i<n;++i)
    if (u[i]>blo) {
      int b = (u[i]>bhi) ? nb : (int)floor(iwid*(u[i]-blo)); // bin # for u[i]
      inv[i] = tempoffs[b]++;              // append to index list in bin b
    }
  for (int i=0;i<n;++i)                    // invert the inv permutation
    if (inv[i]>=0)
      inds[inv[i]] = i;
  
  if (verb) printf("binsort done: %d of %d values to %d+1 bins\n",offs[nb+1],n,nb);
  if (verb>1) {
    printf("inds:\n");
    for (int i=0;i<50;++i)
      cout << inds[i] << endl;
    printf("offs:\n");
    for (int b=0;b<50;++b)
      cout << offs[b] << endl;
  }
    
  return 0;
}
  
int perc3d(int N, double *u, int *d, int* siz, int* nd, int nh, double *hran,
	   double *h0, int verb)
/* 
   Select percolation threshold on gridded 3D real data, from a list of
   thresholds h. Use variant of Newman-Ziff algorithm.

   Inputs: u, a real NxNxN contiguous array of doubles.
           N, the cube side length (positive integer).
	   nh, number of h bins between hran[0] to hran[1].
	   hran, 2-element vector containing lowest and highest h.
	   verb, verbosity (0,1,2,3).
   Outputs:
           d, an int NxNxN domain number array (as in nodal3dziff), at threshold
           siz, a nd-length list of sizes of the domains, at threshold
           nd, the number of domains = max(d), at threshold.
	   h0, estimate of percolation threshold to nearest step.

   Barnett 8/31/17
 */
{
  int n=N*N*N;    // # grid pts
  if (verb) printf("perc3d: N=%d, n=%d, h: [%g,%g] %d levels\n",N,n,hran[0],hran[1],nh);
  double totwid = hran[1]-hran[0];
  double wid = totwid/nh, iwid = 1.0/wid; // h range, bin width & reciprocal
  
  // bin sort...
  vector<int> offs(nh+2);               // bin index starts, last is tot counts
  vector<int> inds(n);                  // bin-ordered list of site indices
  int ier = binsort(n, u, nh, hran[0], iwid, inds, offs, verb);
  int totcts = offs[nh];                // total counts across all bins
  
  vector<int> ptr(n,-1);                // cluster tree parent array
  // join up entire min z face, and max z face
  for (int y=0; y<N; ++y)
    for (int x=0; x<N; ++x) {
      qunion(&ptr[0],0,x+N*y);          // z=0 join to site 0
      qunion(&ptr[0],n-1,x+N*y+n-N*N);  // z=N-1 join to size n-1
    }

  int b=nh;                             // index of bin (start at infinte one)
  int perc = 0;
  while (!perc && b>=0) {               // loop bins in descending order...
    *h0 = hran[0] + b*wid;              // lower end of this bin
    if (verb>2) printf("bin b=%d, h lower=%g\n",b,*h0);
    for (int j=offs[b];j<offs[b+1];++j) {  // all sites in this bin
      int i=inds[j];                    // new site to connect
      if (verb>1)                       // debug
	if (u[i]<*h0 || (b<nh && u[i]>*h0+wid)) printf("%d fail: u=%g\n",i,u[i]);
      int z = (int)(i/(N*N));           // get x,y,z coords of i
      int ixy = i-z*N*N;
      int y = (int)(ixy/N);
      int x = ixy-y*N;
      // use connectedness given by if u>=h, without periodic wrapping:
      if (x>0) {
	int j = i-1;         // x neighbor
	qunion(&ptr[0],i,j);
      }
      if (x<N-1) {
	int j = i+1;         // x neighbor
	qunion(&ptr[0],i,j);
      }
      if (y>0) {
	int j = i-N;         // y neighbor
	qunion(&ptr[0],i,j);
      }
      if (y<N-1) {
	int j = i+N;         // y neighbor
	qunion(&ptr[0],i,j);
      }
      if (z>0) {
	int j = i-N*N;         // z neighbor
	qunion(&ptr[0],i,j);
      }
      if (z<N-1) {
	int j = i+N*N;         // z neighbor
	qunion(&ptr[0],i,j);
      }
    }
    perc = (findroot(&ptr[0],0)==findroot(&ptr[0],n-1));  // faces connected?
    --b;
  }
  if (!perc)
    *h0 = NAN;  // never perced, else *h0 is lower end of bin that first perced
  if (verb) printf("perc done (perc=%d)\n",perc);
  
  // output domain info at threshold, or at lowest h if never perced
  ier = labeldoms(n, &ptr[0], u, d, siz, nd, +1, *h0, verb);  // sets d, siz, nd
  return ier;
}
