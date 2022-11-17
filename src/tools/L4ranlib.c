
#include "../L4common.h"

#ifndef __L4_RAN_LIB__
#define __L4_RAN_LIB__



void L4RandomGeneratorReseed(L4VAR *pv)
{
	time_t RndTimeNow;
	time(&RndTimeNow);
	gsl_rng_set(pv->tkr.gsl_rand_num_gen, 
	    (unsigned long int) (RndTimeNow*(pv->usr.thread_index + 5)));
}


double randnum(L4VAR *pv)
{
	if (pv->tkr.gsl_rand_num_gen == NULL) {
	  const gsl_rng_type * T = gsl_rng_default;
		gsl_rng_env_setup();
		pv->tkr.gsl_rand_num_gen = gsl_rng_alloc (T);
		fprintf(stderr, "L4RandomInitialization %d %lx\n",pv->usr.thread_index, (unsigned long int) pv->tkr.gsl_rand_num_gen);
    pv->tkr.gsl_rng_counter = -999;
  }
  if (pv->tkr.gsl_rng_counter < 0) { // overflow integer
    pv->tkr.gsl_rng_counter = 1;
    L4RandomGeneratorReseed(pv);
  }
  pv->tkr.gsl_rng_counter++;
  const gsl_rng * x = pv->tkr.gsl_rand_num_gen;
  return (double) gsl_rng_uniform_pos( x );
}

double logrand(L4VAR *pv)
{
  return -log(randnum(pv));
}
#endif
