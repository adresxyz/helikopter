#ifndef GPCREGLER
#define GPCREGLER

#define DUMP(varname) printf("%s = \n", #varname);



void CalcRefModelOutput(float y,float w,float W0[]);
void CalcClearArx(ARX *model, float h[]);
void CalcArixOutput(ARX* model,float* y,float* u,float y0[]);
float CalculateGPC(ARX *model, float*y, float *u, float w);
void FillQ(Matrix * mat,float h[]);
void CalcArixOutput(ARX* model,float* y,float* u,float y0[]);
void CalcQT(Matrix* Q,float qT[]);
float CalcArxOutput(ARX* model,float* y,float* u);
void pushStop(float New,float* A,int size);
void pushStart(float New,float* A,int size);
#endif // GPCREGLER

