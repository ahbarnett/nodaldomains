#ifndef DOMAINLIBH
#define DOMAINLIBH

int nodal3dziff(int N, double *u, int *d, int *siz, int *nd, int verb);
int perc3d(int N, double *u, int *d, int* siz, int* nd, int nh, double *h, double *h0, int verb);

#endif
