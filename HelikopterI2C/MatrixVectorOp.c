
#include<math.h>
#include<stdlib.h>
#include<assert.h>
#include <stdio.h>

#include "matrixvectorop.h"
#include "gpcregler.h"
#define DEBUGINFO


Matrix getMatrix(int row,int col, float A)
{
    int i,j;
    float **data = (float **) malloc((unsigned)row * sizeof(float*));
    for( i = 0; i < row; i++) data[i] = (float *)malloc((unsigned) col * sizeof(float));

    for ( i = 0; i < row; ++i) {
        for ( j = 0; j < col; ++j) {
            data[i][j]=A;
        }
    }

    Matrix mat;
    mat.mat = data;
    mat.col = col;
    mat.row = row;
    return mat;
}
Matrix getMatrixEye(int row,int col,float A)
{
    int i,j;
    float ** y = (float **)malloc(row * sizeof(float*));
    for( i = 0; i < row; i++) y[i] = (float *)malloc(col * sizeof(float));

    for ( i = 0; i < row; ++i) {
        for ( j = 0; j < col; ++j) {
            y[i][j]=0;
        }
    }

    for ( i = 0; i < row && i< col; ++i) {
        y[i][i]=A;
    }


    Matrix mat;
    mat.mat = y;
    mat.col = col;
    mat.row = row;
    return mat;
}

void freeMatrix(Matrix* mat)
{
    int i;
    for( i = 0; i < mat->row; i++) free(mat->mat[i]);
    free(mat->mat);
}

void freeMatrixV(Matrix mat)
{
    int i;
    for( i = 0; i < mat.row; i++) free(mat.mat[i]);
    free(mat.mat);
}

void printMatrix(Matrix* mat)
{
    int i,j;
    for(i=0;i<mat->row;++i)
    {
        for(j=0;j<mat->col;++j)
            printf("\t%f",mat->mat[i][j]);
        printf("\n");
    }
    printf("\n");
}

Matrix MatrixMult(Matrix* mat,Matrix* mau)
{
    /// Allocating new Matrix
    Matrix rMatrix;
    /// Copying data
    rMatrix.mat =(float**) matMult(mat->mat,mau->mat,mat->row,mat->col,mau->col);
    rMatrix.col = mau->col;
    rMatrix.row = mat->row;
    return rMatrix;
}

/**
 * @brief MatrixMult
 * @param mat - first matrix
 * @param mau - second matrix
 * @param n - #of rows of f.m.
 * @param m - #of cols of f.m.
 * @param o - #of cols of s.m.
 * @return
 */
float **matMult(float** mat,float** mau,int l,int m,int n )
{

#ifdef DEBUGINFO
    assert(mat!=0);
    assert(mau!=0);
#endif

    int i,j,k;
    float** NewMatrix = (float**)malloc(l*sizeof(float));

    // Allocation of output
    for ( i = 0; i < l; ++i) {
        NewMatrix[i] = (float *)malloc(n * sizeof(float));
    }


    for ( i = 0; i < l; ++i) {
        for ( j = 0; j < n; ++j) {
            NewMatrix[i][j]=0;
            for ( k = 0; k < m; ++k) {
                NewMatrix[i][j]+=mat[i][k] * mau[k][j];
            }
        }
    }

    return NewMatrix;

}

void printTable(float A[],int size)
{
    int i;
    printf("\t");
    for (i = 0; i < size; ++i) {
        printf("%f ",A[i]);
    }
    printf("\n");
}


void moveTableLeft(float* A,int size)
{
    int i;
    for (i = 1; i < size; ++i) {
        A[i-1]=A[i];
    }
}

void pushStop(float New,float* A,int size)
{
    int i;
    for (i = 1; i < size; ++i) {
        A[i-1]=A[i];
    }
    A[size-1]=New;
}

void pushStart(float New,float* A,int size)
{
    int i;
    for (i = size-1; i >0; --i) {
        A[i]=A[i-1];
    }
    A[0]=New;
}
