#include "matrixvectorop.h"
#include "gpcregler.h"
#include "lumatrixinv.h"

#include<stdlib.h>
#include<stdio.h>

#undef SHOW_INFO


#define H 5
#define L 3
#define Alpha 0.3
#define Rho 0.5

///
/// \brief CalculateGPC
/// \param model - model
/// \param yk - obecne wyjście
/// \param w - obecne wymuszenie
/// \return
///
float CalculateGPC(ARX *model, float*y, float *u, float w)
{
    static float w0[H];
    static float h[H];
    static float y0[H];
    static float qT[H];

#ifdef SHOW_INFO
    DUMP(y);printTable(y,model->nA);
    DUMP(u);printTable(u,model->nB+model->k);
#endif


    CalcRefModelOutput(y[model->nA-1],w,w0);

#ifdef SHOW_INFO
    DUMP(w0);printTable(w0,H);
#endif



    CalcClearArx(model,h);
#ifdef SHOW_INFO
    DUMP(h);printTable(h,H);
#endif


    Matrix Q = getMatrix(H,L,0.0);
    FillQ(&Q,h);
#ifdef SHOW_INFO
    DUMP(Q);printMatrix(&Q);
#endif


    CalcQT(&Q,qT);
#ifdef SHOW_INFO
    DUMP(qT);printTable(qT,H);
#endif
    freeMatrix(&Q);

    CalcArixOutput(model,y,u,y0);
#ifdef SHOW_INFO
    DUMP(y0);printTable(y0,H);
#endif

    // Calculating u
    float du = 0;
    int i;
    for ( i = 0; i < H; ++i) {
        du+= qT[i] * (w0[i] - y0[i]);
    }

#ifdef SHOW_INFO
    printf("\n du=%f \n",du);
    printf("\n u=%f \n",du +u[model->nB+model->k-1]);
#endif


    return du +u[model->nB+model->k-1];
}

void CalcQT(Matrix* Q,float qT[])
{
    int i,j,k;
    Matrix QToInv = getMatrix(L,L,0);

    // CalcMatrixToInv
    for ( i = 0; i < L; ++i) {
        for ( j = 0; j < L; ++j) {
            //Inner product
            for (k = 0; k < H; ++k) {
                QToInv.mat[i][j] += Q->mat[k][j] * Q->mat[k][i];
            }
        }
    }

    // Add rho to main diagonal
    for (i = 0; i < L; ++i) {

        QToInv.mat[i][i] += Rho;
    }


#ifdef SHOW_INFO
    DUMP(QToInv);printMatrix(&QToInv);
#endif


    inverse(QToInv.mat,L);


#ifdef SHOW_INFO
    DUMP(QToInv);printMatrix(&QToInv);
#endif

    // Mult QToInv * Q^T = q^T
    for ( k = 0; k < H; ++k) {
        qT[k]=0.0;
        for ( i = 0; i < L; ++i) {
            qT[k]+= QToInv.mat[0][i] * Q->mat[k][i];
        }
    }

    freeMatrix(&QToInv);
}

void FillQ(Matrix * mat,float h[])
{
    int i,j,k;
    for( j=0;j<L;++j)
    {   k=0;
        for( i=j;i<H;++i,++k)
        {
            mat->mat[i][j]=h[k];
        }
    }
}

/**
 * @brief CalcRefModelOutput
 * Obliczenie H-parametrów modelu odniesienia
 * @param y
 * @param w
 * @param W0
 */
void CalcRefModelOutput(float y,float w,float W0[])
{   int i;
    W0[0] = (1.0 - Alpha)*w + Alpha *y;
    for ( i = 1; i < H; ++i) {
        W0[i]=(1.0 - Alpha)*w + Alpha *W0[i-1];
    }
}


/**
 * @brief CalcClearArx
 *
 * Wyliczenie odpowiedzi skokowej z zerowymi warunkami poczatkowymi
 *
 * @param model
 * @param h
 */
void CalcClearArx(ARX* model,float h[])
{
    int j,var=0;
    for (; var < H; ++var) {
        h[var]=0;
        for (j = 0; var-j>0 && j<model->nA; ++j) {
            h[var] -= model->A[j]*h[var-j-1];
        }

        for (j = 0; var - j >=0&& j<model->nB; ++j) {
            h[var] += model->B[j];
        }
    }
}


void CalcArixOutput(ARX* model,float* y,float* u,float y0[])
{
    int i,j;
    float ynext;
    float* wY = (float*) malloc(model->nA *sizeof(float));
    float* wU = (float*) malloc( (model->nB + model->k) *sizeof(float));
    float ui=u[model->nB+model->k-1];


    for(j=0;j<model->nA;j++) wY[j]=y[j];
    for(j=0;j<(model->nB +model->k);j++) wU[j]=u[j];

    moveTableLeft(wU,model->nB+model->k);
    wU[model->nB+model->k-1]=ui;

    // There are H elements of y0 to compute
    for (i = 0; i < H; ++i) {
        ynext=0.0;



        for (j = 0; j<model->nA; ++j) {
            ynext -= model->A[j]*wY[model->nA-j-1];
        }

        for (j = 0;j<model->nB; ++j) {
            ynext += model->B[j]*wU[model->nB-j+model->k-1];
        }

        // move numbers left
        moveTableLeft(wY,model->nA);
        wY[model->nA-1]=ynext;

        moveTableLeft(wU,model->nB+model->k);
        wU[model->nB+model->k-1]=ui;
#ifdef SHOW_INFO
        DUMP(wY);printTable(wY,model->nA);
        DUMP(wU);printTable(wU,model->nB +model->k);
#endif
        y0[i]=ynext;
    }





    free(wU);
    free(wY);

}

float CalcArxOutput(ARX* model,float* y,float* u)
{
    int j;
    float ynext;

    ynext=0.0;

    for (j = 0; j<model->nA; ++j) {
        ynext -= model->A[j]*y[model->nA-j-1];
    }

    for (j = 0;j<model->nB; ++j) {
        ynext += model->B[j]*u[model->nB-j+model->k-1];
    }

    return ynext;
}
