

#ifndef __L4_VectorInitializeFlag
#define __L4_VectorInitializeFlag 1

#define VECTOLLERANCE 1e-8

// typedef struct { int    x, y, z; } I3V;
// typedef struct { double x, y, z; } D3V;
// typedef struct { D3V    x, y, z; } MAT;   // each D3V is a row
// 
I3V IVECZERO = { 0, 0, 0 };
D3V  VECZERO = { 0, 0, 0 };
D3V  VEC_UX  = { 1, 0, 0 };
D3V  VEC_UY  = { 0, 1, 0 };
D3V  VEC_UZ  = { 0, 0, 1 };
MAT  MAT_UNITARY = { 1,0,0, 0,1,0, 0,0,1 };
// 
// #define  PLUS        9001
// #define  MINUS       9002
// #define  TIMES       9003
// #define  CROSS       9005
// #define  CROSSCROSS  9006

double vecdotvec(D3V a, D3V b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
double vecmagsq(D3V a)         { return a.x*a.x + a.y*a.y + a.z*a.z; }
double vecmag(D3V a)           { return sqrt(vecmagsq(a));           }
D3V    vecmuc(D3V a, double b) { a.x *= b; a.y *= b; a.z *= b; return a; }
D3V    vecadc(D3V a, double b) { a.x += b; a.y += b; a.z += b; return a; }

D3V    vecinit(double x, double y, double z) { D3V c = { x, y, z }; return c; } // 7/11/06, 12/26/08
I3V   ivecinit(int x, int y, int z) { I3V c;  c.x = x; c.y = y; c.z = z; return c; }

int    veciszero(D3V a) { return (a.x==0 && a.y==0 && a.z==0 ? 1 : 0); }
int   iveciszero(I3V a) { return (a.x==0 && a.y==0 && a.z==0 ? 1 : 0); }

int    veccompare(D3V a, D3V b) {
  return (fabs(a.x - b.x) > 1e-8 || fabs(a.y - b.y) > 1e-8 || fabs(a.z - b.z) > 1e-8) ? 0 : 1;
}
int   iveccompare(I3V a, I3V b) { return (a.x == b.x && a.y == b.y && a.z == b.z) ? 1 : 0; }

D3V   i3vtod3v(I3V a) { D3V b; b.x = a.x; b.y = a.y; b.z = a.z; return b; }
I3V   d3vtoi3v(D3V a) { I3V b; b.x = a.x; b.y = a.y; b.z = a.z; return b; }
I3V   vecint(D3V a) { I3V b; b.x = a.x; b.y = a.y; b.z = a.z; return b; }
I3V   vecrint(D3V a) { I3V b; b.x = lrint(a.x); b.y = lrint(a.y); b.z = lrint(a.z); return b; }

D3V   vec2unit(D3V a) {
  double mag = vecmag(a);
  return (mag > 0 ? vecmuc(a, 1/mag) : VECZERO);
}  

D3V vecscalar(D3V a, int op, double b)
{ 
  D3V c;
  switch(op) {
    case PLUS:
      c.x = a.x + b; c.y = a.y + b; c.z = a.z + b; break;
    case MINUS:
      c.x = a.x - b; c.y = a.y - b; c.z = a.z - b; break;      
    case TIMES:
      c.x = a.x * b; c.y = a.y * b; c.z = a.z * b; break;
    default: L4stop;
  }
  return c;
}

D3V ivecscalar(I3V a, int op, double b)
{ 
  D3V c;
  switch(op) {
    case PLUS:
      c.x = a.x + b; c.y = a.y + b; c.z = a.z + b; break;
    case MINUS:
      c.x = a.x - b; c.y = a.y - b; c.z = a.z - b; break;      
    case TIMES:
      c.x = a.x * b; c.y = a.y * b; c.z = a.z * b; break;
    default: L4stop;
  }
  return c;
}


D3V vecvec(D3V a, int op, D3V b)
{
  D3V    c;
  switch(op) {
    case PLUS:
      a.x += b.x; a.y += b.y; a.z += b.z; return a;
    case MINUS:
      a.x -= b.x; a.y -= b.y; a.z -= b.z; return a;
    case TIMES:
      a.x *= b.x; a.y *= b.y; a.z *= b.z; return a;
    case CROSS:
      c.x = a.y*b.z - a.z*b.y;
      c.y = a.z*b.x - a.x*b.z;
      c.z = a.x*b.y - a.y*b.x;
      return c;
    case CROSSCROSS:
      return vecvec(a, CROSS, vecvec(a, CROSS, b));
    default: L4stop;
  }
  return VECZERO;
}
I3V ivecvec(I3V a, int op, I3V b)
{
  I3V c;
  switch(op) {
    case PLUS:  c.x = a.x + b.x; c.y = a.y + b.y; c.z = a.z + b.z; break;
    case MINUS: c.x = a.x - b.x; c.y = a.y - b.y; c.z = a.z - b.z; break;      
    case TIMES: c.x = a.x * b.x; c.y = a.y * b.y; c.z = a.z * b.z; break;      
    default: L4stop;
  }
  return c;
}
I3V ivecivec(I3V a, int op, I3V b) { return ivecvec(a, op, b); }



MAT Matinit(double xx, double xy, double xz,  double yx, double yy, double yz,  double zx, double zy, double zz)
{
  MAT m = { xx,xy,xz, yx,yy,yz, zx,zy,zz }; return m;
}

MAT Matini3vec( D3V a, D3V b, D3V c)
{
  MAT m  = { a.x, a.y, a.z,  b.x, b.y, b.z,  c.x, c.y, c.z }; return m;
}

double Matdet(MAT m) {
double mxx = m.y.y*m.z.z - m.y.z*m.z.y;
double mxy = m.y.z*m.z.x - m.y.x*m.z.z;
double mxz = m.y.x*m.z.y - m.y.y*m.z.x;
return m.x.x*mxx + m.x.y*mxy + m.x.z*mxz;
}

MAT Matmuc(MAT m, double c)
{
  MAT n;
  n.x = vecmuc(m.x,c);
  n.y = vecmuc(m.y,c);
  n.z = vecmuc(m.z,c);
  return n;
}

D3V Matxvec(MAT m, D3V v)
{
  D3V c;
  c.x = vecdotvec(m.x,v); c.y = vecdotvec(m.y,v); c.z = vecdotvec(m.z,v);
  return c;
}

MAT Matxmat(MAT a, MAT b)
{
  MAT t = Mattranspose(b);
  MAT c = { vecdotvec(a.x , t.x), vecdotvec(a.x , t.y), vecdotvec(a.x , t.z),
            vecdotvec(a.y , t.x), vecdotvec(a.y , t.y), vecdotvec(a.y , t.z),
            vecdotvec(a.z , t.x), vecdotvec(a.z , t.y), vecdotvec(a.z , t.z) };
  return c;
}

MAT MatPlusMat(MAT a, MAT b)
{
  MAT c = { vecvec(a.x, PLUS, b.x),
            vecvec(a.y, PLUS, b.y),
            vecvec(a.z, PLUS, b.z) };
  return c;
}

MAT MatMinusMat(MAT a, MAT b)
{
  MAT c = { vecvec(a.x, MINUS, b.x),
            vecvec(a.y, MINUS, b.y),
            vecvec(a.z, MINUS, b.z) };
  return c;
}

MAT MatNegative(MAT a)
{
  MAT c = { -a.x.x, -a.x.y, -a.x.z,
            -a.y.x, -a.y.y, -a.y.z,
            -a.z.x, -a.z.y, -a.z.z };
  return c;
}

MAT Mattranspose(MAT m)
{
  MAT t = {  m.x.x,  m.y.x,  m.z.x,
             m.x.y,  m.y.y,  m.z.y,
             m.x.z,  m.y.z,  m.z.z };
  return t;
}

// xx xy xz
// yx yy yz
// zx zy zz
MAT Matinv(MAT m)
{
  MAT i;
  double invdet = 1./Matdet(m);
  i.x.x = (m.y.y*m.z.z - m.y.z*m.z.y)*invdet;   //    xx  yx  zx
  i.x.y = (m.x.z*m.z.y - m.x.y*m.z.z)*invdet;   //    xy  yy  zy
  i.x.z = (m.x.y*m.y.z - m.x.z*m.y.y)*invdet;   //    xz  yz  zz 
  i.y.x = (m.y.z*m.z.x - m.y.x*m.z.z)*invdet;
  i.y.y = (m.x.x*m.z.z - m.x.z*m.z.x)*invdet;
  i.y.z = (m.x.z*m.y.x - m.x.x*m.y.z)*invdet;
  i.z.x = (m.y.x*m.z.y - m.y.y*m.z.x)*invdet;
  i.z.y = (m.x.y*m.z.x - m.x.x*m.z.y)*invdet;
  i.z.z = (m.x.x*m.y.y - m.x.y*m.y.x)*invdet;
  return i;
}


// create unit vectors orthoganal to r (m.x and m.y) and aligned with r (m.z)
MAT createaxialsystem(D3V r)
{
  MAT m;
  m.z = vec2unit(r);
  double absrx = fabs(r.x);
  double absry = fabs(r.y);
  double absrz = fabs(r.z);
  if (absrz < absrx && absrz < absry) {
    m.y = vecvec(m.z, CROSS , VEC_UZ);
  } else if (absrx < absry) {
    m.y = vecvec(m.z, CROSS , VEC_UX);
  } else {
    m.y = vecvec(m.z, CROSS , VEC_UY);
  }
  m.y = vec2unit(m.y);
  m.x = vecvec(m.y, CROSS, m.z);
  return m;
}

//-------------------------------------------------------------------------

#endif
  
  
  
