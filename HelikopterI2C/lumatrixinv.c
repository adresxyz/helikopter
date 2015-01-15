#include "lumatrixinv.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>


void inverse(float **mat, int dim)
{
    int i, j, *indx;
    float **y, d, *col;

    //y = matrix(0, dim - 1, 0, dim - 1);

    y = (float **) malloc(dim * sizeof(float*));
    for( i = 0; i < dim; i++) y[i] = (float *)malloc(dim * sizeof(float));

    for ( i = 0; i < dim; ++i) {
        for ( j = 0; j < dim; ++j) {
            y[i][j]=0;
        }
    }

    indx = (int *)malloc((unsigned)(dim*sizeof(int)));

    //col = vector(0, dim - 1);
    col = (float *)malloc((unsigned)(dim*sizeof(float)));

    ludcmp(mat, dim, indx, &d);
    for (j = 0; j<dim; j++)
    {
        for (i = 0; i<dim; i++) col[i] = 0.0;
        col[j] = 1.0;
        lubksb(mat, dim, indx, col);
        for (i = 0; i<dim; i++) y[i][j] = col[i];
    }
    for (i = 0; i<dim; i++)
        for (j = 0; j<dim; j++)
            mat[i][j] = y[i][j];

    //free_matrix(y, 0, dim - 1, 0, dim - 1);
    // free_vector(col, 0, dim - 1);

    for( i = 0; i < dim; i++){ free(y[i]);}
    free(y);
    free(col);
    free(indx);
}


void ludcmp(float **a, int n, int *indx, float *d)
{
    int i, imax, j, k;
    float   big, dum, sum, temp;
    float   *vv;

    vv = (float*)malloc((unsigned)(n*sizeof(float)));
    if (!vv)
    {
        fprintf(stderr, "Error Allocating Vector Memory\n");
        exit(1);
    }
    *d = 1.0;
    for (i = 0; i<n; i++)
    {
        big = 0.0;
        for (j = 0; j<n; j++)
        {
            if ((temp = fabs(a[i][j])) > big) big = temp;
        }
        if (big == 0.0)
        {
            fprintf(stderr, "Singular Matrix in Routine LUDCMP\n");
            for (j = 0; j<n; j++) printf(" %f ", a[i][j]); printf("/n");
            exit(1);
        }
        vv[i] = 1.0 / big;
    }
    for (j = 0; j<n; j++)
    {
        for (i = 0; i<j; i++)
        {
            sum = a[i][j];
            for (k = 0; k<i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i<n; i++)
        {
            sum = a[i][j];
            for (k = 0; k<j; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }
        if (j != imax)
        {
            for (k = 0; k<n; k++)
            {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = TINY;
        if (j != n - 1)
        {
            dum = 1.0 / a[j][j];
            for (i = j + 1; i<n; i++) a[i][j] *= dum;
        }
    }
    free(vv);
}


void lubksb(float **a, int n, int *indx, float *b)
{
    int i, ip, j, ii = -1;
    float   sum;

    for (i = 0; i<n; i++)
    {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii >= 0)
            for (j = ii; j<i; j++) sum -= a[i][j] * b[j];
        else if (sum) ii = i;
        b[i] = sum;
    }
    for (i = n - 1; i >= 0; i--)
    {
        sum = b[i];
        for (j = i + 1; j<n; j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}
