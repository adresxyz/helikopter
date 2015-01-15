#ifndef MATRIXVECTOROP
#define MATRIXVECTOROP

typedef struct
{
    float** mat;
    int row;
    int col;
}Matrix;

typedef struct
{
    float* A;
    float* B;
    int nA;
    int nB;
    int k;
}ARX;

Matrix MatrixMult(Matrix* mat,Matrix* mau);
Matrix getMatrix(int row,int col, float A);
Matrix getMatrixEye(int row, int col, float A);

float** matMult(float** mat,float** mau,int l,int m,int n );
void freeMatrix(Matrix* mat);

void printMatrix(Matrix* mat);
//void printModel(


void printTable(float[], int size);

void moveTableLeft(float* A,int size);
void pushStop(float New,float* A,int size);
void pushStart(float New,float* A,int size);
#endif // MATRIXVECTOROP

