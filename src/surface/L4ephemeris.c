
//  Data files geocentric_ephemeris_sun_.txt and geocentric_ephemeris_moon.txt
//   are derived from JPL DE405 planetary and lunar ephemerides 
//    downloaded as tables from astropixels.com and concatenated.
//  msg verified by LRO visualization of subsolar point at
//    svs.gsfc.nasa.gov/cgi-bin/details.cgi?aid=4404


double D2R = (M_PI/180);

typedef struct {
  double a0,d0,w0;
} ADW;

typedef struct {
  double phase;
  D3V    emi, mei, msi;  // inertial coord
  D3V    msg, meg;
  double a,d,w;
} L4INERTIALVEC;

// start of ephemeris relative to J2000
// files start 1 jan 2011

//  in makefile
// DEFINITIONS = \
//   "-D YEAR_FIRST=2011"
//   "-D YEAR_LAST=2015"
//   "-D J2000_EPHEMERIS_OFFSET=((2011 - 2000)*365 + 3 - .5)"
//   "-D NUMEPHEMDAYS=(5*365+1)"

L4INERTIALVEC L4inertial[NUMEPHEMDAYS+1];

//---------------------------------------------------------------------------------

D3V radec2vec(double ra, double dc)
{
  D3V p;
  p.x = cos(dc*D2R)*cos(ra*D2R);
  p.y = cos(dc*D2R)*sin(ra*D2R);
  p.z = sin(dc*D2R);
  return p;
}

void earthsunvec(void)
{
  int i;
  char line[200];
  char fn[100];
  sprintf(fn, "%s/lexs-data/solar/sunvec/horizon_ephemeris__sun.txt", getenv("LEXSHOME"));
  FILE *fp = fopen(fn,"r");
  if (fp == NULL) L4stop;
  fgets(line, 199, fp); fgets(line, 199, fp); fgets(line, 199, fp);
  fgets(line, 199, fp); fgets(line, 199, fp);
  double ra, dc;
  for (i = 0; i <= NUMEPHEMDAYS; i++) {
    fscanf(fp,"%*s %*s %lf %lf", &ra, &dc); 
    L4inertial[i].msi = radec2vec(ra,  dc);    // same as earth to sun
  }
  fclose(fp);
}    

void earthmoonvec(void)
{
  int i;
  char line[200];
  char fn[100];
  sprintf(fn, "%s/lexs-data/solar/sunvec/horizon_ephemeris_moon.txt", getenv("LEXSHOME"));
  FILE *fp = fopen(fn,"r");
  fgets(line, 199, fp); fgets(line, 199, fp); fgets(line, 199, fp);
  fgets(line, 199, fp); fgets(line, 199, fp);
  if (fp == NULL) L4stop;
  double ra, dc;
  for (i = 0; i <= NUMEPHEMDAYS; i++) {
    fscanf(fp,"%*s %*s %lf %lf", &ra, &dc); 
    L4inertial[i].emi = radec2vec(ra,  dc);
    L4inertial[i].mei = vecmuc(L4inertial[i].emi, -1.);
  }
  fclose(fp);
}    

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
//  vernal equinox 2000 03 20 07:30 GMT
//   (http://data.giss.nasa.gov/cgi-bin/ar5/srevents.cgi)

ADW updateElements(double d)  //  d = days  from 2000 Jan 1  1200.
{
  double s[14],c[14];
  double E1  = 125.045 -  0.0529921*d;  s[ 1] = sin(E1 *D2R);  c[ 1] = cos(E1 *D2R);
  double E2  = 250.089 -  0.1059842*d;  s[ 2] = sin(E2 *D2R);  c[ 2] = cos(E2 *D2R);
  double E3  = 260.008 + 13.0120009*d;  s[ 3] = sin(E3 *D2R);  c[ 3] = cos(E3 *D2R);
  double E4  = 176.625 + 13.3407154*d;  s[ 4] = sin(E4 *D2R);  c[ 4] = cos(E4 *D2R);
  double E5  = 357.529 +  0.9856003*d;  s[ 5] = sin(E5 *D2R);  c[ 5] = cos(E5 *D2R);
  double E6  = 311.589 + 26.4057084*d;  s[ 6] = sin(E6 *D2R);  c[ 6] = cos(E6 *D2R);
  double E7  = 134.963 + 13.0649930*d;  s[ 7] = sin(E7 *D2R);  c[ 7] = cos(E7 *D2R);
  double E8  = 276.617 +  0.3287146*d;  s[ 8] = sin(E8 *D2R);  c[ 8] = cos(E8 *D2R);
  double E9  =  34.226 +  1.7484877*d;  s[ 9] = sin(E9 *D2R);  c[ 9] = cos(E9 *D2R);
  double E10 =  15.134 -  0.1589763*d;  s[10] = sin(E10*D2R);  c[10] = cos(E10*D2R);
  double E11 = 119.743 +  0.0036096*d;  s[11] = sin(E11*D2R);  c[11] = cos(E11*D2R);
  double E12 = 239.961 +  0.1643573*d;  s[12] = sin(E12*D2R);  c[12] = cos(E12*D2R);
  double E13 =  25.053 + 12.9590088*d;  s[13] = sin(E13*D2R);  c[13] = cos(E13*D2R);
//
  ADW p;
  p.a0 = 269.9949 + (0.0031/36525)*d;
  p.a0 += ( -3.8787*s[ 1]
           - 0.1204*s[ 2]
           + 0.0700*s[ 3]
           - 0.0172*s[ 4]
           + 0.0072*s[ 6]
           - 0.0052*s[10]
           + 0.0043*s[13]);
  p.d0 =  66.5392 + (0.0130/36525)*d;
  p.d0 += (  1.5419*c[ 1]
           + 0.0239*c[ 2]
           - 0.0278*c[ 3]
           + 0.0068*c[ 4]
           - 0.0029*c[ 6]
           + 0.0009*c[ 7]
           + 0.0008*c[10]
           - 0.0009*c[13]);
  p.w0 =  38.3213 + 13.17635815*d - 1.4e-12*d*d;
  p.w0 += (  3.5610*s[ 1]
           + 0.1208*s[ 2]
           - 0.0642*s[ 3]
           + 0.0158*s[ 4]
           + 0.0252*s[ 5]
           - 0.0066*s[ 6]
           - 0.0047*s[ 7]
           - 0.0046*s[ 8]
           + 0.0028*s[ 9]
           + 0.0052*s[10]
           + 0.0040*s[11]
           + 0.0019*s[12]
           - 0.0044*s[13]);
  return p;
}


static double PhaseOffset = 0;
static double LastPhase   = 0;

void AddInertialData(int id) // id is days from 1 Jan 2011 (0 UT)
{
  double day = id + J2000_EPHEMERIS_OFFSET; // days from J2000 (2000/01/01 12:00)
  ADW adw = updateElements(day);
    // Inertial Coordinates = ICRF  (x -> Aries, y is in equatorial plane, and z is north)
    // rotate X(inertial) from aries to ascending node of lunar rotation
  double a = L4inertial[id].a = (90+adw.a0)*D2R;
  double d = L4inertial[id].d = (90-adw.d0)*D2R;
  double w = L4inertial[id].w =  adw.w0*D2R;
    //MAT i2n = {  ca,  sa,   0,     -sa,  ca,   0,      0,   0,   1 };
    // rotate ccw about node to lunar equator (d0 is declination of the lunar pole)
    //MAT n2l = {   1,   0,   0,       0,  cd,  sd,      0, -sd,  cd };
    // rotate x axis to prime meridian about polar axis
    //MAT l2g = {  cw,  sw,   0,     -sw,  cw,   0,      0,   0,   1 };
    // geo coordinate axes in inertial coordinates
    //L4inertial[id].g2i = Matxmat(l2g, Matxmat(n2l, i2n));
  double ca = cos(a), sa = sin(a);
  double cd = cos(d), sd = sin(d);
  double cw = cos(w), sw = sin(w);
  MAT I2G = { ( ca*cw - sa*cd*sw), ( sa*cw + ca*cd*sw), sd*sw ,
              (-ca*sw - sa*cd*cw), (-sa*sw + ca*cd*cw), sd*cw ,
                     sa*sd       ,       -ca*sd       ,  cd     };
  //MAT G2I = Mattranspose(I2G);
  // msg = moon-sun vector in geo coord, i.e. subsolar point
  L4inertial[id].msg = Matxvec(I2G, L4inertial[id].msi);
  L4inertial[id].meg = Matxvec(I2G, L4inertial[id].mei);
  double subsolarlon = atan2(L4inertial[id].msg.y, L4inertial[id].msg.x);
  double subearthlon = atan2(L4inertial[id].meg.y, L4inertial[id].meg.x);
  double subearthlondeg = subearthlon*180/M_PI;
  if (fabs(subearthlondeg) > 9) L4stop;
  double phase = L4fmod2pi(M_PI + subearthlon - subsolarlon);
  if (phase < LastPhase) PhaseOffset += 2*M_PI;
  L4inertial[id].phase = phase + PhaseOffset;
  LastPhase = phase;
}

//---------------------------------------------------------------------------------

int __initInertialData__ = 0;

void L4initializeInertialData(void)
{
  if (__initInertialData__) return;
  __initInertialData__ = 1;
  earthsunvec();
  earthmoonvec();
  int id;
  for (id = 0; id <= NUMEPHEMDAYS; id++) AddInertialData(id);
  // initialize time limits
  //
  //  synodic month = 29.530589 days
  //  draconic year = 29947974.5563 s = 346.620076 days
  //  4 draconic years = (4 + delta)*47*L4SYNODICMONTH
  //  delta = .004205  ( error ~ .1%)
  //  let t_ephem_epoch = 47 synodic months = (1 - .00105) * 4 draconic years
  //    L4NUM_LUNATIONS_EPHEM = 47
  int j, nmcnt = 0;
  double p1, p2 = L4fmod2pi(L4inertial[0].phase);
  for (j = 1; j < NUMEPHEMDAYS; j++) {
    if (nmcnt >= L4NEWMOONCNT) break;
    p1 = p2;
    p2 = L4fmod2pi(L4inertial[j].phase);
    if (p2 < p1) {
      p1 -= 2*M_PI; // p1 < 0
      L4com.NewMoonEphem[nmcnt] = 86400*(j-1 - p1/(p2 - p1));  
      nmcnt++;
    }
  }
  int ilun0 = 11; //  best fit of longitude in  sim --eph data
  L4com.t_ephem_start = L4com.NewMoonEphem[ilun0];
  int ilun1 = ilun0 + L4NUM_LUNATIONS_EPHEM;
  if (ilun1 >= nmcnt) L4stop;
  L4com.t_ephem_limit = L4com.NewMoonEphem[ilun1];
  if (L4com.t_ephem_limit >= 86400*NUMEPHEMDAYS) L4stop;
  L4com.t_ephem_epoch = L4com.t_ephem_limit - L4com.t_ephem_start;
  L4com.t_lunation_avg = L4com.t_ephem_epoch/L4NUM_LUNATIONS_EPHEM; 
}  

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

D3V vecinterpolate( double t, double *p )
{
  double d0 = t/86400;    // ephemeris index = days from 1 Jan 2011
  if (d0 < 0 || d0 >= NUMEPHEMDAYS) L4stop;
  double fq = d0 - (int)d0;
  double fp = 1 - fq;
  double *q = p + sizeof(L4INERTIALVEC)/8;
  D3V v = { fp**p + fq**q,  fp**(p+1) + fq**(q+1),  fp**(p+2) + fq**(q+2) };
  return v;
}

//---------------------------------------------------------------------------------
//     range of t starts at  2011 Jan 1 00:00:00, ends 5 years later
//---------------------------------------------------------------------------------

int L4time2temindex(double t, double *pf1, double *pf2)
{
  double d0 = t/86400;    // ephemeris index = days from 1 Jan 2011
  if (d0 < 0 || d0 >= NUMEPHEMDAYS) L4stop;
  *pf2 = d0 - (int)d0;
  *pf1 = 1 - *pf2;
  return (int)d0;
}

//.................................................................................

MAT L4i2g(double t)
{
  double f1,f2;
  int ix = L4time2temindex(t, &f1, &f2);
  // a = ra + pi/2
  double a1 = L4inertial[ix].a, a2 = L4inertial[ix+1].a;
  if (fabs(a1 - a2) > M_PI) {
    a2 += (a1 > a2) ? 2*M_PI : -2*M_PI;
    if (fabs(a1 - a2) > M_PI) L4stop;
  }
  double a = f1*a1 + f2*a2;
  double ca = cos(a);
  double sa = sin(a);
  // d = pi/2 - dec
  double d1 = L4inertial[ix].d, d2 = L4inertial[ix+1].d;
  if (fabs(d1 - d2) > M_PI) {
    d2 += (d1 > d2) ? 2*M_PI : -2*M_PI;
    if (fabs(d1 - d2) > M_PI) L4stop;
  }
  double d = f1*d1 + f2*d2;
  double cd = cos(d);
  double sd = sin(d);
  // w = longitude from node to zero-meridian
  double w1 = L4inertial[ix].w, w2 = L4inertial[ix+1].w;
  if (fabs(w1 - w2) > M_PI) {
    w2 += (w1 > w2) ? 2*M_PI : -2*M_PI;
    if (fabs(w1 - w2) > M_PI) L4stop;
  }
  double w = f1*w1 + f2*w2;
  double cw = cos(w);
  double sw = sin(w);
  //
  MAT i2g = { ( ca*cw - sa*cd*sw), ( sa*cw + ca*cd*sw), sd*sw ,
              (-ca*sw - sa*cd*cw), (-sa*sw + ca*cd*cw), sd*cw ,
                     sa*sd       ,       -ca*sd       ,  cd     };
  return i2g;
}

MAT L4g2i(double t)
{
  return Mattranspose(L4i2g(t));
}

D3V L4msi(double t)
{
  return vecinterpolate( t, &L4inertial[(int)(t/86400)].msi.x );
}

D3V L4mei(double t)
{
  return vecinterpolate( t, &L4inertial[(int)(t/86400)].mei.x );
}

D3V L4msg(double t)
{
  return vecinterpolate( t, &L4inertial[(int)(t/86400)].msg.x );
}

D3V L4meg(double t)
{
  return vecinterpolate( t, &L4inertial[(int)(t/86400)].meg.x );
}

double L4solphase(double t)
{
  double f1,f2;
  int ix = L4time2temindex(t, &f1, &f2);
  return f1*L4inertial[ ix ].phase + f2*L4inertial[ix+1].phase;
}

double L4solphasemod(double t)
{
  return L4fmod2pi(L4solphase(t));
}

//---------------------------------------------------------------------------------

double L4deltaphase(double t1, double t2)
{
  double dp = L4solphase(t2) - L4solphase(t1);
  while (dp <= 0) dp += 2*M_PI;
  return dp;
}





///////////////////////////////////////////////////////////////////////////////////

void L4ephemeris_test(char *c)
{
//   int testid = atoi(c);
//   L4VAR var[2];
//   L4VAR *pv = &var[0];
//   int id;
//   for (id = 0; id <= NUMEPHEMDAYS; id += 30) {
//     double t = id*86400;
//     if (t >= L4com.t_ephem_limit) exit(0);
//     L4update_ephemeris_parms(pv, t);
//     if (testid == 0) {  // x-geo transformed to inertial should be close to mei
//       D3V xvec = { 1, 0, 0 };
//       D3V mei = Matxvec(pv->trj.g2i, xvec);
//       double dp = atan2(pv->trj.mei.y, pv->trj.mei.x) - atan2(mei.y, mei.x);
//       if (dp >= M_PI) dp -= 2*M_PI;
//       if (dp <  -M_PI) dp += 2*M_PI;
//       printf("%4d  %8.3f\n", id, dp*180/M_PI);
//     } else if (testid == 1) {  // ra and d of z-geo should agree with ephem.
//       D3V zvec = { 0, 0, 1 };
//       D3V zi = Matxvec(pv->trj.g2i, zvec);
//       double sd = zi.z;
//       double cd = sqrt(1 - sd*sd);
//       double d = atan2(sd,cd);  // declination
//       double dd = d - (M_PI/2 - L4inertial[id].d); if (dd > M_PI) dd -= 2*M_PI;
//       double ca = zi.x/cd;
//       double sa = zi.y/cd;
//       double a = atan2(sa,ca);  // right ascension
//       double da = a - (L4inertial[id].a-M_PI/2); if (da < M_PI) da += 2*M_PI;
//       printf("%4d  %8.5f  %8.5f\n",id, da*180/M_PI, dd*180/M_PI);
//     }
//   }
  exit(0);
}  









