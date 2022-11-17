
#include "../L4common.h"
#include "L4vector.c"
#include "L4matpac.c"
#include "L4ranlib.c"

// miscellaneous stuff

double L4fmod2pi( double z )  // same as L4f2modulo(z, 2*M_PI)
{
  if (z >= 0 && z < 2*M_PI) return z;
  double x = z/(2*M_PI);  // integer + fraction
  long   i = x;
  return 2*M_PI*(x >= 0 ? x - i : 1 + x - i);
}

double L4f2modulo(double z, double modulus)
{
  if (z >= 0 && z < modulus) return z;
  double x = z/modulus;  // integer + fraction
  long   i = x;
  return modulus*(x >= 0 ? x - i : 1 + x - i);
}


int L4i2modulo(int ix, int modulus)
{
  return (ix >= 0) ? ix % modulus : (ix + modulus*(1 - ix/modulus));
}

int L4f2index(double z, int modulus, int roundoff)
{
  double x = z/modulus;  // integer + fraction
  long   i = x;
  double f = (x >= 0 ? x - i : 1 + x - i);
  return ( roundoff ? (int)(f*modulus+.5) : (int)(f*modulus) ) % modulus;
}

int L4rad2index(double radian, int modulus, int roundoff)
{
  double x = radian/(2*M_PI);
  long   i = x;
  double f = (x >= 0 ? x - i : 1 + x - i);
  return ( roundoff ? (int)(f*modulus+.5) : (int)(f*modulus) ) % modulus;
}

// Mike Giles   people.maths.ox.ac.uk
// Oxford-Man Institute of Quantitative Finance, Eagle House, Walton Well Road, Ox- ford OX2 6ED

double erfinv(double x)
{
  if (x < 0 || x > 1) {
    fprintf(stderr, " error in erfinv(x), 0 <= x <= 1\n"); exit(0);
  }
  double p, w = -log((1-x)*(1+x));
  if ( w < 5 ) {
    w = w - 2.5;
    p =  2.81022636e-08;
    p =  3.43273939e-07  + p*w;
    p = -3.52338770e-06  + p*w;
    p = -4.39150654e-06  + p*w;
    p =  0.00021858087   + p*w;
    p = -0.00125372503   + p*w;
    p = -0.00417768164   + p*w;
    p =  0.246640727     + p*w;
    p =  1.50140941      + p*w;
  } else {
    w = sqrt(w) - 3;
    p =  -0.000200214257;
    p =   0.000100950558 + p*w;
    p =   0.00134934322  + p*w;
    p =  -0.00367342844  + p*w;
    p =   0.00573950773  + p*w;
    p =  -0.0076224613   + p*w;
    p =   0.00943887047  + p*w;
    p =   1.00167406     + p*w;
    p =   2.83297682     + p*w;
  }
  return p*x;
}
 
void erfinvtest(void)
{
  double x,y,z;
  for (x = .99; x < 1; x += .00001) {
    y = erfinv(x);
    z = erf(y);
    printf("%14.7f  %10.3e  %14.8f %10.2e\n",x,z-x,y, 1-z/x);
  }
  exit(0);
}





