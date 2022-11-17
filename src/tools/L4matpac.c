/* matrix math  */

#include "L4matpac.h"

/* conventions: ara[ Rows ][ Cols ]  */

/**************************************************************************************
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
***************************************************************************************/

#define matval(R,C)       (*(pmat + R*N + C))

void  matadd(double *pA, double *pB, double *pC, int N)
{ int ix; for (ix = 0; ix < N*N; ix++) *pC++ = *pA++ + *pB++; }
void  araadd(double *pA, double *pB, double *pC, int R, int C)
{ int ix; for (ix = 0; ix < R*C; ix++) *pC++ = *pA++ + *pB++; }

void  matsub(double *pA, double *pB, double *pC, int N)
{ int ix; for (ix = 0; ix < N*N; ix++) *pC++ = *pA++ - *pB++; }
void  arasub(double *pA, double *pB, double *pC, int R, int C)
{ int ix; for (ix = 0; ix < R*C; ix++) *pC++ = *pA++ - *pB++; }

void  matadc(double cst, double *pA, double *pC, int N)
{ int ix; for (ix = 0; ix < N*N; ix++) *pC++ = cst + *pA++; }
void  araadc(double cst, double *pA, double *pC, int R, int C)
{ int ix; for (ix = 0; ix < R*C; ix++) *pC++ = cst + *pA++; }

void  matmuc(double cst, double *pA, double *pC, int N)
{ int ix; for (ix = 0; ix < N*N; ix++) *pC++ = cst * *pA++; }
void  aramuc(double cst, double *pA, double *pC, int R, int C)
{ int ix; for (ix = 0; ix < R*C; ix++) *pC++ = cst * *pA++; }

void  matdiag(double cst, double *pU, int N)
{ int ix; double *pu;
  pu = pU; for (ix = 0; ix < N*N; ix++) *pu++ = 0;
  pu = pU; for (ix = 0; ix < N; ix++) { *pu = cst; pu += N+1; }
}
void  matvector2diag(double *pV, double *pU, int N)
{ int ix; double *pu;
  pu = pU; for (ix = 0; ix < N*N; ix++) *pu++ = 0;
  pu = pU; for (ix = 0; ix < N; ix++) { *pu = *pV++; pu += N+1; }
}

void  matxcol(double *pA, double *pB, double *pD, int N)
{
  int iR,iC;
  double sum,*pa;
  for (iR = 0; iR < N; iR++) {
    pa = pA + iR*N;
    sum = 0;
    for (iC = 0; iC < N; iC++) sum += *(pa+iC)**(pB+iC);
    *pD++ = sum;
}}
void  matmul(double *pA, double *pB, double *pD, int N)
{
  int iR,iC,ix;
  double sum,*pa,*pb;
  for (iR = 0; iR < N; iR++) {
    for (iC = 0; iC < N; iC++) {
      pa = pA + iR*N; pb = pB + iC; sum = 0;
      for (ix = 0; ix < N; ix ++)	{ sum += *pb**pa++; pb += N; }
      *pD++ = sum;
}}}
void  aramul(double *pA, double *pB, double *pD, int rAD, int cArB, int cBD)
{
  int iR,iC,ix;
  double sum,*pa,*pb;
  for (iR = 0; iR < rAD; iR++) {
    for (iC = 0; iC < cBD; iC++) {
      pa = pA + iR*cArB; pb = pB + iC; sum = 0;
      for (ix = 0; ix < cArB; ix ++) { sum += *pb**pa++; pb += cBD; }
      *pD++ = sum;
}}}

void   aratranspose(double *pM, double *pT, int rMcT, int cMrT)
{
  int iR,iC;
  double *pt;
  for (iR = 0; iR < rMcT; iR++) {
    pt = pT + iR;
    for (iC = 0; iC < cMrT; iC++) {
      *pt = *pM++; pt += rMcT;
    }
  }
}


/********************************************************************************/

double *pinv,*pmat;
void zeroentry(int rowadj, int coladj, int rowref, int N)
{
  int ic;
  double *pmr, *pma, *pir, *pia, factor;
  pmr = pmat + rowref*N;
  pma = pmat + rowadj*N;
  pir = pinv + rowref*N;
  pia = pinv + rowadj*N;	
  factor = -*(pma + coladj)/(*(pmr + coladj));
  
  if (isinf(factor) || isnan(factor)) L4stop;
  
  for (ic = 0; ic < N; ic++) {
    *pma++ += factor**pmr++;
    *pia++ += factor**pir++;
  }
  *(pma - N + coladj) = 0;
}

int matinv(double *pM, double *pI, int N)
{
  //void    zeroentry(int rowadj, int coladj, int rowref, int N);
  double *pcolmax,*pm,*pi,*p0,*p1;
  int    *irowxref;
  int     iR,iC,ix,iq;
  double  cfn,cmax,qfactor;
  pinv     = pI;
  pmat     = (double *)malloc(N*N*8);
  pcolmax  = (double *)malloc(N*8);
  irowxref = (int    *)malloc(N*4);
  for (iC = 0; iC < N; iC++) {
    *(irowxref + iC) = iC; /* iC is proxy for row index */
    cmax = 0; pm = pM + iC;
    for (iR = 0; iR < N; iR++) {
      if (fabs(*pm) > cmax) cmax = fabs(*pm); pm += N;
    }
    if (cmax == 0) {
      free(pmat); free(irowxref); free(pcolmax);
      fprintf(stderr, "matpac - matinv abort - col = %d\n",iC);
      return -1;
    }
    *(pcolmax + iC) = 1/cmax;
  }
  for (iR = 0; iR < N; iR++) {
    qfactor = 0; 
    ix = *(irowxref + iR);
    for (iC = 0; iC < N; iC++) {
      cfn = fabs(*(pM + ix*N + iC)* *(pM + *(irowxref + iC)*N + iR)* *(pcolmax + iC));
      if (cfn > qfactor) { qfactor = cfn; iq = iC; }
    }
    if (qfactor == 0) {
      free(pmat); free(irowxref); free(pcolmax);
      fprintf(stderr, "matpac - matinv abort - row = %d\n",iR);
      return -2;
    }
    *(irowxref + iR) = *(irowxref + iq);
    *(irowxref + iq) = ix;
  }
  p1 = pmat; pi = pI;
  for (iR = 0; iR < N; iR++) {
    p0 = pM + N**(irowxref + iR);
    ix = *(irowxref + iR);
    for (iC = 0; iC < N; iC++) {
      *p1++ = *p0++;
      *pi++ = (iC == ix ? 1 : 0);
    }
  }
  for (iC = N-1; iC > 0; iC--) {
    for (iR = iC-1; iR >= 0; iR--) {
      if (matval(iR,iC) != 0) zeroentry(iR,iC,iC,N);
    }
  }
  for (iC = 0; iC < N-1; iC++) {
    for (iR = iC+1; iR < N; iR++) {
      if (matval(iR,iC) != 0) zeroentry(iR,iC,iC,N);
    }
  }
  for (iR = 0; iR < N; iR++) {
    cfn = 1./matval(iR,iR);  if (isinf(cfn) || isnan(cfn)) L4stop;
    pi = pI + iR*N;
    for (iC = 0; iC < N; iC++) *pi++ *= cfn;
  }
  free(pmat); free(irowxref); free(pcolmax);
  return 0;
}

/***************  solves Y = MX **************************************************/

int   matsolve(double *pY, double *pM, double *pX, int N)
{
  int ierr;
  double *pI = (double *)malloc(N*N*8);
  if ((ierr = matinv(pM ,pI, N)) == 0) {
    aramul(pI, pY, pX, N,N,1);
  }
  free(pI);
  return ierr;
}

int   arainv(double *pM, double *pI, int rMcT, int cMrT)
{
  int ierr;
  double *pT, *pTM,*piTM;  //,*piTMxT;
  pT     = (double *)malloc(cMrT*rMcT*8);
  pTM    = (double *)malloc(cMrT*cMrT*8);
  piTM   = (double *)malloc(cMrT*cMrT*8);
  //piTMxT = (double *)malloc(cMrT*rMcT*8);
  aratranspose(pM, pT, rMcT, cMrT);
  aramul(pT,pM,pTM, cMrT, rMcT, cMrT);
  if ((ierr = matinv(pTM,piTM, cMrT))) return ierr;
  aramul(piTM  , pT, pI, cMrT, cMrT, rMcT);
  return 0;
}

int   arasolve(double *pY, double *pM, double *pX, int rMrYcT, int cMrXrT, int cYX)
{
  int ierr;
  double *pT,*pTM,*piTM,*piTMxT;
  pT     = (double *)malloc(cMrXrT*rMrYcT*8);
  pTM    = (double *)malloc(cMrXrT*cMrXrT*8);
  piTM   = (double *)malloc(cMrXrT*cMrXrT*8);
  piTMxT = (double *)malloc(cMrXrT*rMrYcT*8);
  aratranspose(pM, pT, rMrYcT, cMrXrT);
  aramul(pT,pM,pTM, cMrXrT, rMrYcT, cMrXrT);
  if ((ierr = matinv(pTM,piTM, cMrXrT)) == 0) {
    aramul(piTM  , pT, piTMxT, cMrXrT, cMrXrT, rMrYcT);
    aramul(piTMxT, pY, pX    , cMrXrT, rMrYcT, cYX   );
  }
  free(pT); free(pTM); free(piTM); free(piTMxT);
  return ierr;
}

/**** W = weight, solves Y = MX where Y=W=[Rx1], M=[RxC], X=[Cx1] *****************/

int  arawsolve(double *pW, double *pY, double *pM, double *pX,
	       int rWYM, int cMrX)
{
  int ierr;
  double *pwm,*pwy,*pww;
  pwm = (double *)malloc(rWYM*cMrX*8);
  pwy = (double *)malloc(rWYM*8);
  pww = (double *)malloc(rWYM*rWYM*8);
  matvector2diag(pW,pww,rWYM);
  aramul(pww,pY,pwy, rWYM,rWYM,1);
  aramul(pww,pM,pwm, rWYM,rWYM,cMrX);
  ierr = arasolve(pwy,pwm,pX,rWYM,cMrX,1);
  free(pwm); free (pwy); free(pww);
  return ierr;
}

/********************************************************************************/
/*------------------------------------------------------------------------
int main(void)
{
  double W[8]    = { 400,400,400,400,400,400,400,400 };
  double Y[8]    = { 1,2,3,4,5,6,7,8 };
  double M[8][4] = { 1,5,3,4,6,9,7,2,  3,2,1,7,6,5,9,8,
                     5,6,2,8,9,3,2,1,  6,7,3,9,2,9,3,2 };
  double X[4];
  arasolve(Y,&M[0][0],X,8,4,1);
  printf(" %8.4f %8.4f %8.4f %8.4f\n",X[0],X[1],X[2],X[3]);
  arawsolve(W,Y,&M[0][0],X,8,4);
  printf(" %8.4f %8.4f %8.4f %8.4f\n",X[0],X[1],X[2],X[3]);

}
-------------------------------------------------------------------------*/
