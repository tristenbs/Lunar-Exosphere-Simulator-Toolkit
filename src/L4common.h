

#ifndef __COMMONINITIALIZED__
#define DEBUG

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/file.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include <gsl/gsl_rng.h>
//#include <gsl/gsl_linalg.h>
#include <pthread.h>

#include "tools/L4vector.h"
#include "tools/L4matpac.h"

#ifdef DEBUG
int __err__, *__errptr__;
#define L4stop   __err__ = (__err__ == 0 ? *__errptr__ : 0) 
#else
#define L4stop  exit(0)
#endif

///////////////////////////////////////////////////////////////////////////////

#define  MAX_THREAD_COUNT   16  // > sysctl  machdep.cpu.thread_count

#define  L4Rm               1737.4
#define  L4GMsun            1.32715e11  // km^3/sec^2
#define  L4GMearth          3.98604e5
#define  L4GMmoon           4902.80031  // GRAIL 
#define  L4E2M              384400.0
#define  L4S2E              1.496e8

#define  L4OMEGA_SUN       (2*M_PI/(365.25*86400))
#define  L4OMEGASQ_SUN     (L4OMEGA_SUN*L4OMEGA_SUN)

#define  L4SYNODICMONTH    (29.530589*86400)
#define  L4SIDEREALMONTH   (27.321661*86400)
#define  L4OMEGA_SYNODIC   (2*M_PI/L4SYNODICMONTH)
#define  L4OMEGA_SIDEREAL  (2*M_PI/L4SIDEREALMONTH)
#define  L4OMEGASQ_MOON    (L4OMEGA_SIDEREAL*L4OMEGA_SIDEREAL)
#define  L4DRACONICYEAR     29947974.5563

#define  AEON              (3.16e7*1e9)

#define  L4KB               1.38e-26     // 1e-10*erg/deg  (vel in km/s)
#define  L4PLANCKSCONST     6.6251e-27   // erg sec
#define  L4KOVERH          (1e10*L4KB/L4PLANCKSCONST)
#define  L4HOVERK          (L4PLANCKSCONST/(1e10*L4KB))
#define  L4ADSTEPS_PER_LUN  180
#define  L4DTADS           (L4SYNODICMONTH/L4ADSTEPS_PER_LUN)

#define  L4N_DSRB_ITER       12

// ephemeris wraparound = 4 draconic years ~ 47 lunations
#define  L4NUM_LUNATIONS_EPHEM   47


#define solarconstant        1.361e3      // W / m^2
#define stefanboltzmannconst 5.670373e-8  // W / m^2 / K^−4
#define albedo               0.13         // DIVINER working number
#define emissivity_diviner   0.98         // DIVINER working number

#define T4SUBSOLAR ( (solarconstant/stefanboltzmannconst)*(1 - albedo)/emissivity_diviner )
#define T4RAD                3.2e5        // 1.8e-6 W cm^-2 = 23.88^4/sigma

/* Mass is relative to AMU2GRAMS */
#define  AMU2GRAMS        1.660538921e-24
#define  MELE             (9.1083e-28/AMU2GRAMS)
// http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=html&isotype=some
#define  HYDROGEN           1.00782503207
#define  HELIUM             4.00260325415
#define  CARBON            12.0000000
#define  NITROGEN          14.0030740048
#define  OXYGEN            15.99491461956
#define  NEON              19.9924401754
#define  ARGON36           35.967545106
#define  ARGON40           39.9623831225
#define  WATER            ( 2*HYDROGEN + OXYGEN )
#define  HYDROXYL         (   HYDROGEN + OXYGEN )
#define  HYDROGEN2        ( 2*HYDROGEN )
#define  METHANE          (4*HYDROGEN + CARBON)
#define  CARBONMONOXIDE   (CARBON +   OXYGEN)
#define  CARBONDIOXIDE    (CARBON + 2*OXYGEN)
#define  SILICON           28.0855
#define  TITANIUM          47.867
#define  ALUMINIUM         26.9815386
#define  FE_54             53.9396105
#define  FE_56             55.9349375
#define  FE_57             56.9353940
#define  MAGNESIUM         24.3050
#define  CALCIUM           40.078 
#define  SODIUM            22.98976928
#define  MANGANESE         54.938045
#define  POTASSIUM         39.0983

#define PHOTOION_NEON      8.52048e6
#define PHOTOION_ARGON     1.6e6
#define PHOTOION_CARBON    (1/.41e-6)
#define PHOTOION_NITROGEN  (1/.19e-6)
#define PHOTOION_OXYGEN    (1/.21e-6)

// from Huebner et al. 1992 (notebook 30 Oct 2012)
// quiet sun
#define PHOTOION_QS_H_  (1 / 7.3e-8)
#define PHOTOION_QS_He  (1 / 5.2e-8)
#define PHOTOION_QS_C_  (1 / 4.1e-7)
#define PHOTOION_QS_O_  (1 / 2.1e-7)
#define PHOTOION_QS_Ar  (1 / 3.1e-7)
// active sun
#define PHOTOION_AS_H_  (1 / 1.7e-8)
#define PHOTOION_AS_He  (1 / 1.5e-7)
#define PHOTOION_AS_C_  (1 / 9.2e-7)
#define PHOTOION_AS_O_  (1 / 5.9e-7)
#define PHOTOION_AS_Ar  (1 / 6.9e-7)

#define ATOMH_      0
#define ATOMHE      1

enum { ATOMO_, ATOMSI, ATOMTI, ATOMAL, ATOM54, ATOM56,
       ATOM57, ATOMMG, ATOMCA, ATOMNA, ATOMMN, ATOMK_,
       NATOM } enumGrainAtomIndex;
typedef struct {
  char   name[8];
  double m, Z;
} ATOM;

enum {
  _BULKA, _BULK, _A11, _A12, _A14, _A15, _A16, _A17, _L46, _L40, _L44, _FAN, NSOILDEF
} soil_composition_list;

typedef struct { int y, d, h, m; } YDHM;
typedef struct {
  char *shortlabel, *datalabel, *swlabel;
  YDHM t1;
  YDHM t2;
} SWFILEDEF;

#define LADEEPROXY  0
#define APOLLO17    1

#define NUMSWSS    20


// n = 6e23 atom/mole / 23 g/mole * 2.6 g/cc
// sigma = pi*1e-8^2
// free path = 1/n*sigma = 4.7A
// atom sep. = n^-1/3 = 2.45A
// depth of reflection = 1/2 atom sep. ~1/4 free path
// used in L4scatter.c and L4solarwindscatter.c
//#define  ELASTICSCATTERREFLECTOR  -3  // -0.25 


#define  GrainDensity      2.7

#define  L4ALT_MAX         11.0
#define  L4SURFACE_MAX    (L4Rm + L4ALT_MAX)
#define  L4ALT_MIN         -9.1
#define  L4SURFACE_MIN    (L4Rm + L4ALT_MIN)

#define  L4SURF_UNINITIALIZED      -1
#define  L4SUNLIT                   1
#define  L4SHADOW                   0

#define  L4REGOLITHEXITNOP         0x300
#define  L4REGOLITHEXIT            0x301

enum {
  REG_IMPACT, REG_RADIOGENIC, REG_DIFFUSION, REG_TEMSPEC 
          //  , REG_POPPED // sputter release is negligible
} mode_inregolith;

#define L4ision   (pv->trj.qion  == 0 ? 0 : 1)
#define L4wasion  (pv->trx.qion  == 0 ? 0 : 1)

// bit flags for L4trj.state
#define  L4_SW_PRI                  1
#define  L4_SW_SEC                  2
#define  L4_EXO                     4
#define  L4_RADICAL                 8

// event designator 
#define  L4NEWATOM                  0
#define  L4PHASEERROR               1
#define  L4SW_ESCAPE                2
#define  L4SW_TRAP                  3
#define  L4EX_ESCAPE                4
#define  L4EX_TRAP                  5
#define  L4UNIMPLANT                6
#define  L4TIMETRAP                 7
#define  L4GARDENUNTRAP             8
#define  L4PHOTOLYSIS               9
#define  L4SOLWINDCEX              10
#define  L4PHOTOSCATT              11
#define  L4NUMEVENT                12

#define  L4SW_NOTRACK              20
#define  L4IMPLANT                 21
#define  L4FORCEDTRAP              25
#define  L4SW_NO_IMPACT            26
#define  L4MSG_TRAPPED             0x201
#define  L4MSG_FROCED_ABORT        0x201
#define  L4MSG_FORCED_UNTRAP       0x202
#define  L4MSG_FORCED_TRAP         0x204
#define  L4MSG_FORCED_EXIT         0x205

#define  L4FORCED_TRACKER_ABORT    0x999c

// SourceType for L4SurfaceSource()
#define L4TRK_TEMSPEC             111
#define L4TRK_DIFFUSION           112
#define L4TRK_RADIOGENIC          113
#define L4TRK_EXOSPHERE           114
#define L4TRK_SOLARWIND           115
#define L4TRK_CHARGE_EXCHANGE     116

// options for trj.bdry
#define BDRY_NOP         0  //  not in cell or xing out2in or out2out
#define SOLAR_EVT        1
#define BDRY_AUDIT       2
#define BDRY_WRAP        4
#define BDRY_CELL        8
#define BDRY_IS_IN      32
#define BDRY_WAS_IN     64

#define TIME_GARDEN_UNTRAP       (86400*365.25*1e8)
#define TIME_GARDEN_UNIMPLANT    (86400*365.25*2e8)
#define TIME_COLD_TRAP           (86400*365.25*1e9)

typedef void (*L4AUDITEV)( void * , int, int );
typedef void (*L4AUDITSV)( void * , int, ... );
typedef void (*L4AUDITFN)( void * );
typedef void (*L4INT2INT)( void * );

#define  L4NEWMOONCNT            62
#define  DEFAULT_ROUGHNESS_CFN  0.1

typedef struct {
  int             atom_def;
  int             catalog_index;
  // LADEE era
  int             YearFirst;  // L4dissociation.c
  int             DayFirst ;  // L4dissociation.c
  int             YearLast ;  // L4dissociation.c
  int             DayLast  ;  // L4dissociation.c
  int             SoilComposition;
  int             GeotailHWDeg;
  int             T_Debye;
  int             vel_dist_order;
  int             diffusion_pop;
  int             verbose;          // not implemented
  L4AUDITEV       EventAuditor;
  L4AUDITFN       SurfaceImpactAuditor;
  L4AUDITFN       SurfaceExitAuditor;
  L4AUDITFN       ShellAuditor;
  L4AUDITFN       CellAuditor;
  L4AUDITSV       UpdateCellBoundaryIndicies;
  // primary parameters copied from atom def
  double          MassAMU;
  int             Z;
  int             tem_total_accom;
  double          dt_max;
  double          RadEsc;
  double          PhotonScatterTime;
  double          PhotonWavelengthNm;
  double          SWTrapProbability;
  double          IonInelasticityCfn;
  double          IonSurvivalProb;
  double          ElasticScatterReflector;
  // default parameters 
  double          InterstitialBarrierEnergy;
  double          ImpactUntrapTem;
  double          infrared_cfn;
  double          R;
  double          Cal_Per_Mole_Tref;
  double          Tref;
  double          nu_A;
  double          nu_C;
  double          nu_D;
  double          nu_q;
  double          nu_x1;
  double          MaxIonEV;
  double          T_low_limit;
  double          log_t_garden;
  double          nu_ydsrb_lim;
  int             T_water_trap;
  int             nop;
  // calculated
  double          t_lunation_avg;
  double          t_ephem_start;
  double          t_ephem_limit;
  double          t_ephem_epoch;
  double          InterstitialBarrierUmagSq;
  double          InterstitialBarrierUmag;
  double          GeotailBoundaryLo;
  double          GeotailBoundaryHi;
  double          PhotonScatterDeltaV;
  double          U_desorb;
  double          accommodation_cfn;
  double          t_garden;
  double          NewMoonEphem[L4NEWMOONCNT];
} L4COMDEF;

typedef struct{
  double   MassAMU;
  int      Z;
  int      nop1;
  // inertial coordinates
  D3V      Gravity;
  D3V      r;                           // radius vector (km)  +x is anti-solar, +z is north pole
  D3V      u;                           // velocity vector (km/s)
  //
  double   phase;                       // phase of moon (new moon = 0 phase)
  double   radius;                      // moon-centered radius
  double   radsq;                       // radius-squared
  double   ur;                          // u DOT r
  double   usq;                         // u^2
  double   umag;                        // |u|
  double   tof;                         // time of ballistic, zeroed at boundary crossings
  double   t_ephem;
  //
  double   energyproxy;                 //  v^sq - GM/r
  //
  int      state;
  int      qion;
  int      superescape;
  unsigned int bdry;
  int      insw;
  int      insun;
} L4TRAJECTORY;

// sub-particle surface point -- independent of phase
typedef struct {
  MAT      geo2reg, reg2geo;
  D3V      n;     // fractal normal
  D3V      r;
  D3V      u;
  double   n_dot_ruv;
  double   radius;
  double   theta;
  double   sintha;
  double   costha;
  double   phi;
  double   sinphi;
  double   cosphi;
  int      valid;
  int      intrap;
} L4GEO_DATA;

// subsurface
typedef struct {
  MAT    reg2grn, grn2reg;
  D3V    ureg;        // intra-grain vel; surface normal is z axis
  D3V    ugrn;        // grain surface normal is z axis
  double umag;        // inter/intra-grain speed
  double zreg;        // depth of penetration in regolith
  double zgrn;			  // depth of penetration in grain
  double T;
  double vibvelsq;
  double t_encounter;
} L4REG_DATA;


// shadow and insolation parameters
typedef struct {
  double ifr_cfn;
  double sol_cfn;
} L4SOL_DATA;
  
typedef struct {
  double    LogDevPhotolysis;
  double    LogDevPhotonScatter;
  double    LogDevSW_encounter;
  double    tads;
	const gsl_rng *gsl_rand_num_gen;
	int       gsl_rng_counter;
	int       jday;
} L4TRACKER;

typedef struct {
  char            name[8];
  double          MassAMU;
  int             Z;
  int             tem_total_accom;
  // exosphere
  double          PhotonScatterTime;
  double          PhotonWavelengthNm;
  double          CexCrossSection;
  double          SWTrapProbability;
  double          IonInelasticityCfn;
  double          IonSurvivalProb;
  double          ElasticScatterReflector;
} L4ATOMDEF;


#define SWSTACKSIZE  50
typedef struct {
    L4TRAJECTORY trj;
    L4GEO_DATA   geo;
    L4SOL_DATA   sol;
    L4REG_DATA   reg;
  	L4TRACKER    tkr;
} IMPLANTSTACK;

typedef struct {
  IMPLANTSTACK  ImplantLocator;
  D3V           UntrapUgrn[SWSTACKSIZE];
  double        Untrapzgrn[SWSTACKSIZE];
  int           UntrapCount;
  int           Implant;
} L4STACKS;

// user defined data --------------------------------

typedef struct {
	int            thread_index;
	int            create;
	int            catalog_index;
  int            itsrc;
  int            itim, xtim;
  int            ilon, xlon;
  int            ilat, xlat;
  int            irad, xrad;
} L4USR;
//---------------------------------------------------


typedef struct {
  L4USR          usr;
  L4USR          usx;
  L4TRAJECTORY   trj, trx;
  L4GEO_DATA     geo;
  L4SOL_DATA     sol;
  L4REG_DATA     reg;
  L4TRACKER      tkr;
  L4STACKS       stk;
} L4VAR;

//////////////////////////////////////////////////////////////////////////////////
#ifdef __L4RUN_MODULE
  void   EventAuditorNOP( void *pv, int state, int event ) { }
  void   SurfaceImpactAuditorNOP( void *pv) { }
  void   SurfaceExitAuditorNOP( void *pv) { }
  void   ShellAuditorNOP( void *pv ) { }
  void   CellAuditorNOP( void *pv ) { }

  SWFILEDEF L4SW_Properties[NUMSWSS] = {
    "L\0" , "LDE\0" , "_LADEE.bin\0"    ,  2013,335, 0, 0,   2014, 60, 0, 0,
    "l\0" , "LAD\0" , "_LADEE_HCO.bin\0",  2014, 30, 0, 0,   2014,100, 0, 0,
    "a\0" , "A17\0" , "_Apollo17.bin\0" ,  1973,  1, 0, 0,   1973,267, 0, 0
  };


  L4COMDEF  L4com = {
                             4, // int             atom_def;
                             0, // int             catalog_index;
                                // LADEE era
                          2013, // int             YearFirst;
                           330, // int             DayFirst;
                          2014, // int             YearLast;
                           105, // int             DayLast;
                         _BULK, // int             SoilComposition;
        									  30, // int             GeotailHWDeg;
                           450, // int             T_Debye;
                             3, // int             vel_dist_order;
                             0, // int             diffusion_pop;
                             0, // int             verbose;
      			  &EventAuditorNOP, // L4AUDITEV       EventAuditor;
      &SurfaceImpactAuditorNOP, // L4AUDITFN       SurfaceImpactAuditor;
        &SurfaceExitAuditorNOP, // L4AUDITFN       SurfaceExitAuditor;
              &ShellAuditorNOP, // L4AUDITFN       ShellAuditor;
               &CellAuditorNOP, // L4AUDITFN       CellAuditor;
                          NULL, // L4AUDITSV       UpdateCellBoundaryIndicies;
                                // primary parameters copied from atom def
                       ARGON40, // double          MassAMU;
                             2, // int             Z;
                             0, // int             tem_total_accom;
                          2000, // double          dt_max;
                       20*L4Rm, // double          RadEsc; // approx. neutral potential (L1,L2)
                             0, // double          PhotonScatterTime;
                             0, // double          PhotonWavelengthNm;
                            .1, // double          SWTrapProbability;
                            .5, // double          IonInelasticityCfn;
                           .09, // double          IonSurvivalProb;
                          -.25, // double          ElasticScatterReflector;
                                // default parameters 
                             3, // double          InterstitialBarrierEnergy;
  					      				 500, // double          ImpactUntrapTem;
  									     0.001, // double          infrared_cfn; // cf. ir.c
                     1.9872041, // double          R;            // Cal/K/mole
                          7300, // double          Cal_Per_Mole_Tref;
                           100, // double          Tref;
                          .919, // double          nu_A;
                          .465, // double          nu_C;
                          .610, // double          nu_D
                          .500, // double          nu_q;
                           .85, // double          nu_x1;
                           250, // double          MaxIonEV;
                             0, // double          T_low_limit;
                           6.0, // double          log_t_garden; // 6Ma
                           4.0, // double          nu_ydsrb_lim;
                             0, // int             T_water_trap;
                             0  // int             nop;
  };
  ATOM Atom[NATOM] = {
    "Oxy  \0",  OXYGEN    ,   8,
    "Si   \0",  SILICON   ,  14,
    "Ti   \0",  TITANIUM  ,  22,
    "Al   \0",  ALUMINIUM ,  13,
    "Fe-54\0",	FE_54     ,  26,
    "Fe-56\0",  FE_56     ,  26,
    "Fe-57\0",  FE_57     ,  26,
    "Mg   \0",  MAGNESIUM ,  12,
    "Ca   \0",  CALCIUM   ,  20,
    "Na   \0",  SODIUM    ,  11,
    "Mn   \0",  MANGANESE ,  25,
    "K    \0",  POTASSIUM ,  19
  };
  double      CumulativeAbundance[NATOM];
  double      MineralAtomConc[NATOM];
  double      MeanAtomicMass;
#else
  extern SWFILEDEF  L4SW_Properties[NUMSWSS];
  extern L4COMDEF   L4com;
  extern ATOM       Atom[NATOM];
  extern double     CumulativeAbundance[NATOM];
  extern double     MineralAtomConc[NATOM];
  extern double     MeanAtomicMass;

  void   EventAuditorNOP( void *pv, int state, int event );
  void   SurfaceImpactAuditorNOP( void *parm );
  void   SurfaceExitAuditorNOP( void *parm , int);
  void   ShellAuditorNOP( void *pv );
  void   CellAuditorNOP( void *pv );

#endif

//////////////////////////////////////////////////////////////////////////////////

// tracking subroutine called by front-end code
int    L4Tracker(L4VAR *, int Create, int New);
// Surface Physics subroutines called by tracker
void   L4ParameterInitialization(L4ATOMDEF *);
void   L4VibEnergyAvg_1d( L4VAR *pv );
void   L4getSeleneLocation(L4VAR *, int);
void   L4Terminators(L4VAR *pv);
double L4getRadius(double theta, double phigeo);
double ThermalVelocity(L4VAR *, double Tem);
void   ReadRegolithPixelArray();

int    L4InRegolith(L4VAR *, int);
void   L4geo2trj(L4VAR *pv);
void   L4reg2geo(L4VAR *pv);
void   L4trj2geo(L4VAR *pv);
void   L4geo2reg(L4VAR *pv);
void   L4geo_surf_transformation(L4VAR *pv);
void   L4reg2grn(L4VAR *pv);
void   L4grn2reg(L4VAR *pv);

void   L4UpdateGravity(L4VAR *pv);
void   L4UpdateTrajectoryData(L4VAR *pv, double dt);

// utility subroutines
void   L4initSurfaceData(void);
void   L4SoilComposition(void);
int    L4SolarWindIonRandomLocation(L4VAR *);
int    L4SolarWindIonFixedLocation(L4VAR *);
int    L4UpdateThaPhiIndicies(L4VAR *);
void   L4InitializeTrjIx(L4VAR *);
void   L4DiurnalTemperature(L4VAR *, double lat_deg, double lon_deg, double dlon_deg);
double L4GeographicTemperature(L4VAR *, double lat_deg, double lon_deg, double phase_deg);
int    L4AtomGrainEncounter(L4VAR *);
int    L4SW_Encounter(L4VAR *, int retrapflag);
int    L4ElasticScatterBarrierTest(L4VAR *);
int    L4TargetId(L4VAR *);
void   L4DesorbVel(L4VAR *);
D3V    L4oscillator(L4VAR *, double mass);
void   L4PhononInitialization(void);
void   L4init_velocitydist();
double L4debye_oscillator(double mass, double T);
int    L4SolarWindIonVelocity(L4VAR *, double InitialPhase, int PosFlag, ...);
void   L4initSolarWindData(SWFILEDEF *);
void   L4readUntrapBufferStack(L4VAR *, int);
int    L4pushUntrapBufferStack(L4VAR *, D3V, double);
void   L4getImplantLocator(L4VAR *);
void   L4putImplantLocator(L4VAR *);
void   L4clearAllStacks(L4VAR *);
void   L4InitializeRunTimeParameters(void);
void   getFileName(char *nameaftersrc, char *fullname);

double CellTemperature(int ir, int ip);

void   L4InitializeSurfaceParms(L4VAR *);
void   L4trj_primary(L4VAR *);
int    getRandomSolarWindIon(L4VAR *, int, int);
D3V    L4v_cross_Be2m(L4VAR * );

double getMeanSolarWindParms(L4VAR *);
void   L4ChargeExchangeEvent(L4VAR *);
void   L4TrapBalanceTrackerLoop(L4VAR *pv);
void   L4GeoTrapBalanceTrackerLoop(L4VAR *pv, double, double, double);

void   L4oscillatortest(L4VAR *pv);

void   L4RandomInitialization(int ix);
double randnum(L4VAR *);
double logrand(L4VAR *);

int    L4i2modulo(int ix, int modulus);
double L4f2modulo(double z, double modulus);
double L4fmod2pi(double);
int    L4rad2index(double radian, int modulus, int roundoff);
int    L4f2index(double x, int modulus, int roundoff);
double erfinv(double x);

double L4theta2ilat(double tha);

int    LolaSunshineDetector(D3V rgeo, D3V msg);

void   PolarPofile( L4VAR *pv, char *c );

void   init_get_nu(L4ATOMDEF *patm, int i);
void   L4init_nu_suntan_parms(void);

int    L4BallisticStep(L4VAR *pv);

void   L4_InitDissocData(void);
void   L4_Photolysis(L4VAR *pv);

void   L4initializeInertialData(void);
int    L4time2temindex(double t, double *pf1, double *pf2);
void   L4update_ephemeris_parms(L4VAR *pv, double t);
D3V    L4msg(double t);
D3V    L4meg(double t);
double L4msg_dot_n(L4VAR *pv, double t);
double L4msg_dot_r(L4VAR *pv, double t);
double L4solphase(double t);
double L4solphasemod(double t);

void   L4update_solar_data(L4VAR *pv);

void   L4initTemperatureData(void);
double L4Temperature(L4VAR *, int);
void   L4LADEEtemperaturearray(L4VAR *pv);

D3V    L4mei(double t);
D3V    L4msi(double t);
D3V    L4msg(double t);
MAT    L4i2g(double t);
MAT    L4g2i(double t);
double L4phigeo(L4VAR *, double);

void   L4geo2trj_r(L4VAR *pv);
void   L4geo2trj_v(L4VAR *pv);

void   L4ephemeris_test(char *c);
void   L4Tconvolution_test(void);
double L4PhotolysisTime(L4VAR *pv);

void   EphemEpochTest(char *c);

double Ar_SW_encounter_deviate(L4VAR *pv);
int    Ar_SW_scatter_result(L4VAR *pv);
int    inSolarWindTest(L4VAR *pv);
void   L4PhotonScatterEvent(L4VAR *pv);
double L4get_nu(L4VAR *pv);

//////////////////////////////////////////////////////////////////////////////////
#define __COMMONINITIALIZED__
#endif
