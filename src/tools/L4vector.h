
#ifndef __L4_VECTOR_INIT
#define __L4_VECTOR_INIT

typedef struct { short int x, y, z; } S3V;
typedef struct { int       x, y, z; } I3V;
typedef struct { double    x, y, z; } D3V;
typedef struct { D3V       x, y, z; } MAT;   // each D3V is a row

extern I3V IVECZERO;
extern D3V  VECZERO;
extern D3V  VEC_UX ;
extern D3V  VEC_UY ;
extern D3V  VEC_UZ ;
extern MAT  MAT_UNITARY;

#define  PLUS        9001
#define  MINUS       9002
#define  TIMES       9003
#define  CROSS       9005
#define  CROSSCROSS  9006

double vecdotvec(D3V a, D3V b);
double vecmagsq(D3V a);
double vecmag(D3V a);
D3V    vec2unit(D3V a);

D3V    vecmuc(D3V a, double b);
D3V    vecadc(D3V a, double b);

D3V    vecinit(double x, double y, double z);
I3V   ivecinit(int x, int y, int z);

int    veciszero(D3V a);
int   iveciszero(I3V a);

int    veccompare(D3V a, D3V b);
int   iveccompare(I3V a, I3V b);

D3V   i3vtod3v(I3V a);
I3V   d3vtoi3v(D3V a);
I3V   vecint(D3V a);
I3V   vecrint(D3V a);

D3V vecscalar(D3V a, int op, double b);
D3V ivecscalar(I3V a, int op, double b);
D3V vecvec(D3V a, int op, D3V b);
I3V ivecvec(I3V a, int op, I3V b);
I3V ivecivec(I3V a, int op, I3V b);
MAT Matinit(double xx, double xy, double xz,  double yx, double yy, double yz,  double zx, double zy, double zz);
MAT Matini3vec( D3V a, D3V b, D3V c);
double Matdet(MAT m);
MAT Matmuc(MAT m, double c);
D3V Matxvec(MAT m, D3V v);
MAT Matxmat(MAT a, MAT b);
MAT MatPlusMat(MAT a, MAT b);
MAT MatMinusMat(MAT a, MAT b);
MAT MatNegative(MAT a);
MAT Mattranspose(MAT m);
MAT Matinv(MAT m);
MAT createaxialsystem(D3V r);

//  from matpac.c

int   matinv(double *pM, double *pI, int N);                                         // I = M_inv
int   arasolve(double *pY, double *pM, double *pX, int rMrYcT, int cMrXrT, int cYX); // X = M^-1 x Y
void  aratranspose(double *pM, double *pT, int rMcT, int cMrT);                      // T = M_transpose
int   arawsolve(double *pW, double *pY, double *pM, double *pX, int rWYM, int cMrX);

void  matadd(double *pA, double *pB, double *pC, int N);
void  matsub(double *pA, double *pB, double *pC, int N);
void  matmul(double *pA, double *pB, double *pC, int N);                             // C = A x B
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
