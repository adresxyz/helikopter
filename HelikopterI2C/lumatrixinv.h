#ifndef LUMATRIXINV
#define LUMATRIXINV
#define   TINY 1.0e-20

void inverse(float **mat, int dim);
void ludcmp(float **a, int n, int *indx, float *d);
void lubksb(float **a, int n, int *indx, float *b);
#endif // LUMATRIXINV

