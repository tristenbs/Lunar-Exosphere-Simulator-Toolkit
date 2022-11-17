
#ifndef __MATPAC_HEADER_FILE__
#define __MATPAC_HEADER_FILE__

int   matinv(double *pM, double *pI, int N);
int   arasolve(double *pD, double *pA, double *pC, int rDrA, int cArC, int cDC);
void  aratranspose(double *pM, double *pT, int rMcT, int cMrT);
int   arawsolve(double *pW, double *pY, double *pM, double *pX, int rWYM, int cMrX);

void  matadd(double *pA, double *pB, double *pC, int N);
void  matsub(double *pA, double *pB, double *pC, int N);
void  matxcol(double *pA, double *pB, double *pD, int N);
void  matmul(double *pA, double *pB, double *pC, int N);
void  matadc(double cst, double *pA, double *pC, int N);
void  matmuc(double cst, double *pA, double *pC, int N);
void  araadd(double *pA, double *pB, double *pC, int R, int C);
void  arasub(double *pA, double *pB, double *pC, int R, int C);
void  aramul(double *pA, double *pB, double *pD, int rAD, int cArB, int cBD);
void  araadc(double cst, double *pA, double *pC, int R, int C);
void  aramuc(double cst, double *pA, double *pC, int R, int C);

void  matdiag(double cst, double *pU, int N);
void  matvector2diag(double *pV, double *pD, int N);

#endif
