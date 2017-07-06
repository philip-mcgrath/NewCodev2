#ifndef TRAJECTORY
#define TRAJECTORY

#include <gsl/gsl_rng.h>

double ddd4, ddd, ppower, de; //indexed by trajectory index
double **RR, **PP, *dhat, *f, *dgam, *Pperp /*perp element of momentum to d - position dependent*/;
double abs_d /*trajectory dependent*/, Pdotdhat /*component of momentum parallel to d (function)*/;

double sina, cosa, alpha; //trajectory dependent

double  abszsum0, argzsum0;
double *abszsum1, *argzsum1, **realsum, **imagsum; //0 - summand trajectory dependent ; 1 - sum
double  habszsum0, hargzsum0;
double *habszsum1, *hargzsum1, **hrealsum, **himagsum;


#endif
