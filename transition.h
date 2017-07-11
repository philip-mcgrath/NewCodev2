#ifndef TRANSITION
#define TRANSITION

#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   <vector>

#include "functions.h"

void setwww();
double dens_init_0(double *x,double *p, int k);
double dens_init_1(double *x,double *p, int k);
double dens_init_2(double *x,double *p, int k);
double dens_init_3(double *x,double *p, int k);

double obs_0(double *x,double *p, int k);
double obs_1(double *x,double *p, int k);
double obs_2(double *x,double *p, int k);
double obs_3(double *x,double *p, int k);

double H_0(double *x,double *p, int k);
double H_1(double *x,double *p, int k);
double H_2(double *x,double *p, int k);
double H_3(double *x,double *p, int k);


#endif
