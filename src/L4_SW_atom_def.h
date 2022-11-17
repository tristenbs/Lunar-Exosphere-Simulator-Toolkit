
#ifndef __SWCOMMONINITIALIZED__
#define __SWCOMMONINITIALIZED__ 1


#define U_HYDROGEN 0
#define U_HELIUM   1
#define U_NEON     2
#define U_AR36     3
#define U_ARGON    4
#define NATOMDEF   8

#ifdef __L4SW_MODULE

L4ATOMDEF L4sw_atom_def[NATOMDEF] = {
//  [0]
                 "__h\0", // char            name[8];
                HYDROGEN, // double          MassAMU;
                       1, // int             Z
                     500, // int             tem_total_accom;
                1/1.3e-3, // double          PhotonScatterTime;
                  121.57, // double          PhotonWavelengthNm;
                   2e-15, // double          CexCrossSection; // calculated
                     .10, // double          SWTrapProbability;
                     .50, // double          IonInelasticityCfn;
                     .09, // double          IonSurvivalProb;
                    -.25, // double          ElasticScatterReflector;
//  [1]
                 "_he\0", // char            name[8];
                  HELIUM, // double          MassAMU;
                       2, // int             Z
                     500, // int             tem_total_accom;
                   2.9e4, // double          PhotonScatterTime;
                    58.4, // double          PhotonWavelengthNm;
                 2.7e-17, // double          CexCrossSection; // Stier and Barnett, 1956
                     .06, // double          SWTrapProbability;
                     .50, // double          IonInelasticityCfn;
                     .09, // double          IonSurvivalProb;
                      -3, // double          ElasticScatterReflector;
//  [2]
                 "_ne\0", // char            name[8];
                    NEON, // double          MassAMU;
                      10, // int             Z
                     500, // int             tem_total_accom;
                       0, // double          PhotonScatterTime;
                       0, // double          PhotonWavelengthNm;
                 2.2e-16, // double          CexCrossSection;
                     .10, // double          SWTrapProbability;
                     .50, // double          IonInelasticityCfn;
                     .09, // double          IonSurvivalProb;
                      -3, // double          ElasticScatterReflector;
//  [3]
                 "_ar\0", // char            name[8];
                 ARGON36, // double          MassAMU;
                      18, // int             Z
                     500, // int             tem_total_accom;
                       0, // double          PhotonScatterTime;
                       0, // double          PhotonWavelengthNm;
                1.03e-15, // double          CexCrossSection; Phelps et al. 1992
                     .10, // double          SWTrapProbability;
                     .50, // double          IonInelasticityCfn;
                     .09, // double          IonSurvivalProb;
                      -3, // double          ElasticScatterReflector;
//  [4]
                 "_ar\0", // char            name[8];
                 ARGON40, // double          MassAMU;
                      18, // int             Z
                     500, // int             tem_total_accom;
                       0, // double          PhotonScatterTime;
                       0, // double          PhotonWavelengthNm;
                1.03e-15, // double          CexCrossSection; Phelps et al. 1992
                     .10, // double          SWTrapProbability;
                     .50, // double          IonInelasticityCfn;
                     .09, // double          IonSurvivalProb;
                      -3, // double          ElasticScatterReflector;
//  [5]
                 "ch4\0", // char            name[8];
                 METHANE, // double          MassAMU;
                      10, // int             Z
                     500, // int             tem_total_accom;
                       0, // double          PhotonScatterTime;
                       0, // double          PhotonWavelengthNm;
                       0, // double          CexCrossSection;
                       1, // double          SWTrapProbability;
                     .50, // double          IonInelasticityCfn;
                     .09, // double          IonSurvivalProb;
                      -3, // double          ElasticScatterReflector;
//  [6]
                "_co\0", // char            name[8];
          CARBONMONOXIDE, // double          MassAMU;
                      14, // int             Z
                     500, // int             tem_total_accom;
                       0, // double          PhotonScatterTime;
                       0, // double          PhotonWavelengthNm;
                       0, // double          CexCrossSection;
                       1, // double          SWTrapProbability;
                     .50, // double          IonInelasticityCfn;
                     .09, // double          IonSurvivalProb;
                      -3, // double          ElasticScatterReflector;
//  [7]
                 "co2\0", // char            name[8];
           CARBONDIOXIDE, // double          MassAMU;
                      22, // int             Z
                     500, // int             tem_total_accom;
                       0, // double          PhotonScatterTime;
                       0, // double          PhotonWavelengthNm;
                       0, // double          CexCrossSection;
                       1, // double          SWTrapProbability;
                     .50, // double          IonInelasticityCfn;
                     .09, // double          IonSurvivalProb;
                      -3  // double          ElasticScatterReflector;
};

#else
extern L4ATOMDEF L4sw_atom_def[];
#endif
#endif

