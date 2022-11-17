
/*===========================================================================*/
/* calorie = 4.184 joules */
//#define  THERMALCONDUCTIVITY   (1e-5*1e7)       // W/cm/K * erg/s/W  = erg/cm/sec/K (Keihm and Langseth, 73
//#define  Sigma                  5.679e-5        // Stefan-Boltzman constant,  erg/cm2/s/K^4



//--------------------------------------------------------------------------------

// cf.lexs-data/TemperatureData/simtem.c

#define L4SECPERHOUR  (L4SYNODICMONTH/24)

#define CVORDER     4
#define POWCV      .5   // .25  //.5
#define OFSCV       1   // 1.0  //.5

#define CV_STEP   (360/L4ADSTEPS_PER_LUN)  // 2 degrees per convolution step 
#define T4_SIZE   (180/CV_STEP + 1)
#define CV_SIZE   (3*360/CV_STEP)
#define CV_DT      L4DTADS                 // (CV_STEP*L4SYNODICMONTH/360)

static double  cv_array[CV_SIZE];
static double  T4subsolar = T4SUBSOLAR;
static double  T4lolim;

static int     __initTemperatureDataFlag__ = 0;


//--------------------------------------------------------------------------------

//#define CV_INV_TEST  1

void L4initTemperatureData(void)
{
  if (__initTemperatureDataFlag__) return;
  __initTemperatureDataFlag__ = 1;
  int i,j,n;
  T4lolim = pow(L4com.T_low_limit, 4.0);
// get Langseth model
  double  T, T4Langseth[T4_SIZE];
  char fn[100];
  sprintf(fn, "%s/lexs-data/temperaturedata/temperature_data.txt",
                getenv("LEXSHOME"));
  FILE *fp = fopen(fn,"r");
  if (fp == NULL) L4stop;
  while (fscanf(fp,"%d  %lf", &i, &T) != EOF) {
    if (i > 270) break;
    if (((i-90) % CV_STEP) == 0) {
      j = (i-90)/CV_STEP;
      if (j >= T4_SIZE) L4stop;
      T4Langseth[j] = T*T*T*T - T4lolim;
    }
  }
  fclose(fp);
//
  double G[CV_SIZE][CVORDER];
  for (j = 0; j < CV_SIZE; j++) {
    for (n = 0; n < CVORDER; n++) G[j][n] = pow((j+OFSCV), -(n+1)*POWCV);
  }
//
  double M[T4_SIZE][CVORDER];
  for (i = 0; i < T4_SIZE; i++) {
    for (n = 0; n < CVORDER; n++) M[i][n] = 0;
    for (j = 0; j < CV_SIZE; j++) { 
      double sdotr = sin((j-i+OFSCV)*CV_STEP*M_PI/180);  // time measured back from sunset
      if (sdotr > 0) {
        for (n = 0; n < CVORDER; n++) M[i][n] += sdotr*G[j][n];
      }
    }
  }
  double cv_cfns[CVORDER];
  // int arasolve(double *pY, double *pM, double *pX, int rMrYcT, int cMrXrT, int cYX)
  if (arasolve(&T4Langseth[0], &M[0][0], &cv_cfns[0], T4_SIZE, CVORDER, 1)) L4stop;
  // void  aramul(double *pA, double *pB, double *pD, int rAD, int cArB, int cBD)
  aramul(&G[0][0], &cv_cfns[0], &cv_array[0], CV_SIZE, CVORDER, 1);
  //
  for (j = 0; j < CV_SIZE; j++) {
    if (cv_array[j] < 0) {
      for (j = 0; j < CV_SIZE; j++) printf("%4d  %10.0f\n",j,cv_array[j]);
      exit(0);
    }
  }
  double t4_dark_noon = 0;
  for (j = 0; j < CV_SIZE; j++) {
    double sdotr = cos((j+OFSCV)*CV_STEP*M_PI/180);  // time measured from noon
    if (sdotr > 0) t4_dark_noon += sdotr*cv_array[j];
  }
  t4_dark_noon += T4lolim;
  T4subsolar = T4SUBSOLAR*T4SUBSOLAR/(T4SUBSOLAR + t4_dark_noon + T4RAD);
}

void L4Tconvolution_test(void)
{
  int it;
  for (it = 0; it < 90; it += 10) {
    L4com.T_low_limit = it;
    __initTemperatureDataFlag__ = 0;
    L4initTemperatureData();    
  }
  exit(0);
}

//-----------------------------------------------------------------------------------
double cv_sdotn[MAX_THREAD_COUNT][CV_SIZE];


double L4T4_dark(L4VAR *pv, int adsorb_step)
{
  double *psdn0 = &cv_sdotn[pv->usr.thread_index][0];
  int j;
  if (adsorb_step == 0) {
    double *ps = psdn0;
    double dt  = CV_DT*(OFSCV-1);
    for (j = 0; j < CV_SIZE; j++) {
      dt          += CV_DT;
      D3V    msg   = L4msg(pv->trj.t_ephem - dt);
      double sdotn = vecdotvec(msg, pv->geo.n);
      if (sdotn <= 0) {
        *ps++ = 0;
      } else {
        *ps++ = (LolaSunshineDetector(pv->geo.r, msg)) ? sdotn : 0;
      }
    }
  } else {
    double *ps2 = psdn0 + (CV_SIZE-1);
    double *ps1 = ps2 - 1;
    for (j = (CV_SIZE-1); j > 0; j--) *ps2-- = *ps1--;
    double dt    = CV_DT*OFSCV;
    D3V    msg   = L4msg(pv->trj.t_ephem - dt);
    double sdotn = vecdotvec(msg, pv->geo.n);
    if (sdotn <= 0) {
      *psdn0 = 0;
    } else {
      *psdn0 = (LolaSunshineDetector(pv->geo.r, msg)) ? sdotn : 0;
    }
  }
  double t4_dark = 0;
  for (j = 0; j < CV_SIZE; j++) t4_dark += cv_array[j]**psdn0++;
  double t4sol = t4_dark + T4lolim*pv->geo.sintha;
  if (t4sol < 0) L4stop;
  return t4sol;
}

//--------------------------------------------------------------------------------

double L4Temperature(L4VAR *pv, int adstep)
{
  if (!pv->geo.valid) L4stop;
  double t4sol = L4T4_dark(pv, adstep);
  D3V msg = L4msg(pv->trj.t_ephem);
  double sdotn = vecdotvec(msg, pv->geo.n);
  if (sdotn > 0) {
    if (LolaSunshineDetector(pv->geo.r, msg)) t4sol += T4subsolar*sdotn;
  } else {
    double sdotr = vecdotvec(msg, pv->geo.r);
    if (sdotr > 0) t4sol += pv->sol.ifr_cfn*T4subsolar*sdotr/pv->geo.radius;
  }
  pv->reg.T = pow((t4sol + T4RAD),0.25);
  if (pv->reg.T > 400 || pv->reg.T < 23) L4stop;
  if (isnan(pv->reg.T)) L4stop;
  if (fabs(pv->geo.costha) < .09 && pv->reg.T < 60) L4stop;
  return pv->reg.T;
}

//--------------------------------------------------------------------------------


