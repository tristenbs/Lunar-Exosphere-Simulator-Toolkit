
#define POLE_THETA_DEG     25
#define NPOLE_MAXT      ( POLE_THETA_DEG*64 )
char    pole_maxt[2][2*NPOLE_MAXT][2*NPOLE_MAXT];
double  pole_sintha, pole_dxy;

void initPoleMaxT(void)
{
  float maxt[2*NPOLE_MAXT][2*NPOLE_MAXT];
  int i;
  FILE *fp;
  char fn[100];
  sprintf(fn, "%s/Tmax/maxT-N.bin", getenv("LEXSCURRENT"));
  fp = fopen(fn,"r");
  fread(&maxt[0][0], 4, (2*NPOLE_MAXT*2*NPOLE_MAXT), fp);
  fclose(fp);
  float *pf = &maxt[0][0];
  char  *pc = &pole_maxt[0][0][0];
  for (i = 0; i < (2*NPOLE_MAXT*2*NPOLE_MAXT); i++) {
    *pc++ = (*pf > 0 && *pf <= L4com.T_water_trap) ? 'x' : '_';
    pf++;
  }
  sprintf(fn, "%s/Tmax/maxT-S.bin", getenv("LEXSCURRENT"));
  fp = fopen(fn,"r");
  fread(&maxt[0][0], 4, (2*NPOLE_MAXT*2*NPOLE_MAXT), fp);
  fclose(fp);
  pf = &maxt[0][0];
  pc = &pole_maxt[1][0][0];
  for (i = 0; i < (2*NPOLE_MAXT*2*NPOLE_MAXT); i++) {
    *pc++ = (*pf > 0 && *pf <= L4com.T_water_trap) ? 'x' : '_';
    pf++;
  }
  pole_sintha = sin(POLE_THETA_DEG*M_PI/180);
  pole_dxy = NPOLE_MAXT/pole_sintha;
}

//--------------------------------------------------------------------------------

#define LL_DXDN         (M_PI/(180*64))
#define LL_DNDX         ( 64*180/M_PI)

#define LL_THA           (180*64)
#define LL_PHI           (360*64)

double  LL_sinx[LL_PHI+1], LL_cosx[LL_PHI+1];
float   LL_alt[LL_THA+1][LL_PHI+1];

#define HE_THA            1440
#define HE_PHI            2880
short int hurst_exp[HE_THA][HE_PHI];

//double LL_slope[36][90];   // degrees - average slope, pix = 5 deg lat, 4 deg lon.
int __SurfaceDataInit__ = 0;

//-------------------------------------------------------//
//                                                       //
//      lola_64S.bin is derived from lola_64.bin         //
//         by lola_interp.c                              //
//                                                       //
//-------------------------------------------------------//



void L4initSurfaceData(void)
{
  if (__SurfaceDataInit__) return;
  __SurfaceDataInit__ = 1;
  initPoleMaxT();
  int ip;
  FILE *fp;
  char fn[100];
  sprintf(fn, "%s/lexs-data/lola/lola_64S.bin", getenv("LEXSHOME"));
  fp = fopen(fn, "r");
  fread(&LL_alt[0][0], 4, (LL_THA+1)*(LL_PHI+1) , fp);
  fclose(fp);
  sprintf(fn, "%s/lexs-data/lola/hurst/HE_8.IMG", getenv("LEXSHOME"));
  fp = fopen(fn, "r");
  fread(&hurst_exp[0][0], 2, HE_THA*HE_PHI, fp);
  fclose(fp);
  for (ip = 0; ip <= LL_PHI; ip++) {
    LL_sinx[ip] = sin(ip*LL_DXDN);
    LL_cosx[ip] = cos(ip*LL_DXDN);
  }
  LL_sinx[0] = LL_sinx[LL_PHI] = 0;
  LL_cosx[0] = LL_cosx[LL_PHI] = 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

double hurst_exponent(L4VAR *pv)
{
  int it = L4rad2index(pv->geo.theta , HE_THA, 1);
  int ig = L4rad2index(pv->geo.phi, HE_PHI, 1);
  return hurst_exp[it][ig]*.000025 + .5;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
#define  L4DELTARADMIN (4* (2*M_PI*L4Rm/(360*64))) // 4 x 0.4738 = 1.895 km resolution

int LolaSunshineDetector(D3V rgeo, D3V msg)
{
  // ray-trace to sun in geo-coord
  int    it,ip;
  double radius, ct, p;
  D3V    dr = { msg.x*L4DELTARADMIN, msg.y*L4DELTARADMIN, msg.z*L4DELTARADMIN };
  for ( ; ; ) {
    rgeo.x += dr.x;
    rgeo.y += dr.y;
    rgeo.z += dr.z;
    radius = sqrt(rgeo.x*rgeo.x + rgeo.y*rgeo.y + rgeo.z*rgeo.z);
    if (radius >= L4SURFACE_MAX) break;
    ct = rgeo.z/radius;
    it = acos(ct)*LL_DNDX + .5;
    p  = atan2(rgeo.y, rgeo.x);
    ip = L4rad2index(p, LL_PHI, 1);
    if (radius < (L4Rm + LL_alt[it][ip])) return 0;
  }
  return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

D3V unitvec(int it, int ig)
{
  return vecinit( LL_sinx[it]*LL_cosx[ig], LL_sinx[it]*LL_sinx[ig], LL_cosx[it] );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

    //          NORTH                 SOUTH
    //   nw _|_________|_ ne   nw _|___\|___\|_ ne
    //       |\        |           |\   |\   |
    //       |  \ f=+1 |           | \  | \  |
    //       |    \    |           |  \ |  \ |
    //       | f=-1 \  |       sw _|___\|___\|_ se
    //   sw _|____m___\|_ se   nw  |\   m    |  ne
    //   nw  |\   |\   |  ne       |  \      |
    //       | \  | \  |           |    \    |
    //       |  \ |  \ |           |      \  |
    //   sw _|___\|___\|_ se   sw _|________\|_ se
    //       |    |    |           |    |    |
    //      igw  igm  ige         igw  igm  ige

////////////////////////////////////////////////////////////////////////////////////////////////////////

void L4getSeleneLocation(L4VAR *pv, int mode) 
// mode=0: trj.r and t_ephem,  no geo vars are valid
//      1: geo.theta, costha, .sintha, .phi, .sinphi, .cosphi are valid
//      2: geo.costha and .phi, other geo vars are invalid
{
  if (mode == 0) {           // derive from trj.r, radius
    D3V r = Matxvec(L4i2g(pv->trj.t_ephem), pv->trj.r);
    // inverse test
    double rmag = vecmag(r);
    if (fabs(rmag - pv->trj.radius) > 1e-10) {
      printf("%8.3f %8.3f %8.3f\n", rmag, pv->trj.radius, vecmag(pv->trj.r));
      L4stop;
    }
    pv->geo.costha = r.z/rmag;
    pv->geo.sintha = sqrt(1 - pv->geo.costha*pv->geo.costha);
    pv->geo.theta  = atan2(pv->geo.sintha, pv->geo.costha);
    double rxy     = rmag*pv->geo.sintha;
    pv->geo.cosphi = r.x/rxy;
    pv->geo.sinphi = r.y/rxy;
    pv->geo.phi    = atan2(pv->geo.sinphi, pv->geo.cosphi);
    if (pv->geo.phi < 0) pv->geo.phi += 2*M_PI;
  } else if (mode == 2) {             // only geo.costha and geo.phi are specified
    pv->geo.sintha = sqrt(1 - pv->geo.costha*pv->geo.costha);
    pv->geo.theta  = atan2(pv->geo.sintha, pv->geo.costha);
    pv->geo.cosphi = cos(pv->geo.phi);
    pv->geo.sinphi = sin(pv->geo.phi);
  }
  D3V uvecrgeo = { pv->geo.sintha*pv->geo.cosphi,
                     pv->geo.sintha*pv->geo.sinphi,
                       pv->geo.costha };
  int itn, its, igw, ige;
  D3V rnw, rne, rsw, rse,  unw, use;
  //
  // north-south indicies
  itn = pv->geo.theta*LL_DNDX; if (itn >= LL_THA) {
    if (itn == LL_THA ) itn = LL_THA-1; else L4stop;
  }
  its = itn + 1;
  igw   = pv->geo.phi*LL_DNDX;
  ige   = igw + 1; if (ige > LL_PHI) L4stop;
  rnw   = vecmuc((unw = unitvec(itn, igw)), L4Rm + LL_alt[itn][igw]);
  rne   = vecmuc(       unitvec(itn, ige) , L4Rm + LL_alt[itn][ige]);
  rsw   = vecmuc(       unitvec(its, igw) , L4Rm + LL_alt[its][igw]);
  rse   = vecmuc((use = unitvec(its, ige)), L4Rm + LL_alt[its][ige]);
  // float precision =  7 significant decimal digits i.e. altitude = xx.xxxxx
  int fractal = (vecdotvec(uvecrgeo, vecvec(unw, CROSS, use)) > 0) ? 1 : -1;
  // (unw CROSS use)  points toward north-east
  D3V r0, v1, v2;
  if (fractal == 1) {
    //if (itn ==   0   ) L4stop; else fprintf(stderr, "north pole -----------\n");
    r0   = rse;
    v1   = vecvec(rne, MINUS , r0);
    v2   = vecvec(rnw, MINUS , r0);
  } else {
    //if (its == LL_THA) L4stop; else fprintf(stderr, "south pole------------\n");
    r0   = rnw;
    v1   = vecvec(rsw, MINUS , r0);
    v2   = vecvec(rse, MINUS , r0);
  }
  pv->geo.n = vec2unit(vecvec(v1, CROSS, v2));
  pv->geo.n_dot_ruv = vecdotvec(pv->geo.n, uvecrgeo);
  if (pv->geo.n_dot_ruv < 0) L4stop;
  // equation of plane: Ax + By + Cz = D
  // (Ax+By+Cz) = (nxyz dot rxyz) = radius x (nxyz dot uvecrgeo)
  // D = (nxyz dot r0) = constant
  // rxyz = radius * uvecrgeo ( located on plane )
  //
  pv->geo.radius = vecdotvec(pv->geo.n, r0)/pv->geo.n_dot_ruv;
  pv->geo.r      = vecmuc(uvecrgeo, pv->geo.radius);
  
  //
  if (pv->geo.radius < 1000) L4stop;
  //fprintf(stderr, "sloc  %5d %5d   %5d %5d  %10.4f %8.5f\n", itP,itE, igm,fractal, pv->trj.radius - pv->geo.radius, pv->trj.ur);
  //
  pv->sol.ifr_cfn       = L4com.infrared_cfn;  //*hurst_exponent(pv);
  pv->sol.sol_cfn       = 1;   /// 1 - pv->sol.ifr_cfn;
  pv->geo.valid         = 1;
  if (pv->geo.sintha < pole_sintha) {
    int ix = NPOLE_MAXT + pole_dxy*uvecrgeo.x;
    int iy = NPOLE_MAXT + pole_dxy*uvecrgeo.y;
    int ip = (uvecrgeo.z > 0) ? 0 : 1;
    pv->geo.intrap = (pole_maxt[ip][ix][iy] == 'x') ? 1 : 0;
  } else pv->geo.intrap = 0;
}
//   T4_sun_noon      = T4SUBSOLAR * pv->sol.sun_cfn;
//   T4_ifr_noon      = T4SUBSOLAR * pv->ifr.ifr_cfn;

////////////////////////////////////////////////////////////////////////////////////////

double L4maxalt( double tha1, double tha2, double phi1, double phi2)
{
  int ip1 = (int)(phi1*LL_DNDX);
  int ip2 = (int)(phi2*LL_DNDX);
  int it1 = (int)(tha1*LL_DNDX);
  int it2 = (int)(tha2*LL_DNDX);
  double max = -9999;
  int ip,it;
  for (ip = ip1; ip <= ip2; ip++) {
    int iP = (ip + LL_PHI) % LL_PHI;
    for (it = it1; it <= it2; it++) {
      if (LL_alt[it][iP] > max) max = LL_alt[it][iP];
    }
  }
  return max;
}

/////////////////////////////////////////////////////////////////////////////////////


