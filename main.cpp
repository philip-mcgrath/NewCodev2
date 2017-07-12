#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <complex>


#include "Global.h"
#include "RandomNumberGenerator.h"
#include "MultiPathProcessing.h"



// Global Data
extern PathInfoQueue_t  path_info_queue;
extern PathDataVector_t multi_paths_data;

using namespace std;



char  datafilename[80];
FILE   *stream;
// ========================  INPUT VALUES  ======================== //
int N_bath;
int N_slice;
int Ncut;
double timestep;
double T;
unsigned long init_seed;
int Nsample;
double w_max;
double eta;
int beta;
double delta;
double ppower;

// ============================================================================
// Variables
// ============================================================================
double ddd4;
double ddd;
double Dt;
double Pdotdhat;
double de;
double abs_d;
double sina;
double cosa;



const gsl_rng_type * TT;
gsl_rng * rr;

complex<double> z = 1.0;
complex<double> I(0, 1);

// ============================================================================
// Vectors
// ============================================================================

double *mww;
double *mu;
double *sig;
double *dtdtm;
double *dgam;
double *dhat;
double *R1;
double *v;
double *f;
double *c;
double *m;
double *w;



// Pointers to Functions ====================================
double (*dens_init[4])(double*, double*, int );
double (*obs[4])(double*, double*, int );
double (*obs1[4])(double*, double*, int );

double (* www[6][4][4])();
void (*force[4])(double *);

// =============================================================================
// Multi Path Processing Program
// =============================================================================
// 0 - 1 - 5  - 21
//            - ...
//       - 6
//       - 7
//       - 8
//   - 2 - 9
//       - 10
//       - 11
//       - 12
//   - 3 - 13
//       - 14
//       - 15
//       - 16
//   - 4 - 17
//       - 18
//       - 19
//       - 20
// breath first processing

int main() {

    // ================================================================================================================
    // INPUT
    // ================================================================================================================

    /* Initialize Stream - Scope 0
    cout << "Print information about new stream: "<< endl;
    cout << "Input datafilename, N_bath, N_slice, Ncut" << endl;
    cout << "timestep, T, init_seed, Nsample" << endl;
    cout << "w_max, eta, beta, delta, power" << endl;
    cin >> datafilename >> N_bath >> N_slice >> Ncut >> timestep >> T >> init_seed >> Nsample >> w_max >> eta >> beta >> delta >> ppower;
    */
    /* initialize stream  - scope 0 */
    cout << " Print information about new stream:" << endl;
    cout << "Input datafilename" << endl;
    cin >> datafilename;
    N_bath = 200;
    N_slice = 20;
    Ncut = 10;
    timestep = 0.05;
    T = 15;
    init_seed = 0;
    Nsample = 100;
    w_max = 3;
    eta = 0.13;
    beta = 25;
    delta = 0.8;
    ppower = 100000;


    // ================================================================================================================
    // Memory Allocation
    // ================================================================================================================
    mww = new double[N_bath];
    mu = new double[N_bath];
    sig =  new double[2*N_bath];
    dtdtm = new double[N_bath];
    dgam = new double[N_bath];
    dhat = new double[N_bath];
    R1 = new double[N_bath];
    v = new double[N_bath];
    f = new double[N_bath];
    c = new double[N_bath];
    m = new double[N_bath];
    w = new double[N_bath];



    // Random Number Generator
    // fixed seed for reproducibility, otherwise use RandomState()
    cout << "Root Seed: " << init_seed << endl;
    RandomState random_state = RandomState(init_seed); // root path

    // !!! Multi Paths Data
    long n_data1D   = 4; // dimension parameters // inital values needed by children
    long n_data2D_1 = 4; // sum1 length
    long n_data2D_2 = N_slice; // abszum, argzum, habszum, hargzsum
    long n_paths = pow(N_PATHS, (N_LEVELS+1.0)) - 1;

    multi_paths_data.resize(n_paths, PathData(n_data1D, n_data2D_1, n_data2D_2));

    // ================================================================================================================
    // Initialization Values
    // ================================================================================================================

    gsl_rng_env_setup();

    TT = gsl_rng_default;
    rr = gsl_rng_alloc(TT);


    dens_init[0] = dens_init_0; dens_init[1] = dens_init_1;
    dens_init[2] = dens_init_2; dens_init[3] = dens_init_3;
    obs[0] = obs_0; obs[1] = obs_1; obs[2] = obs_2; obs[3] = obs_3;
    obs1[0] = H_0; obs1[1] = H_1; obs1[2] = H_2; obs1[3] = H_3;
    ddd4 = delta*delta*0.25;
    ddd =  delta*delta;
    Dt  = T/N_slice;


    bath_para(eta,w_max);       /* compute system parameters etc */

   //  bath corresponds to eq. 53//   for (i = 0; i < N_bath; i++)
    for (int i = 0; i < N_bath; ++i)
        mu[i] = beta*w[i]*0.5;
    for (int i = 0; i < N_bath; ++i){
        sig[i] = 1.0/sqrt(w[i]*2.0*tanh(mu[i]));
        mww[i] = -m[i]*w[i]*w[i];
        dtdtm[i] = -0.5*timestep*timestep/m[i];
    }
    for (int i = 0; i < N_bath; ++i)
        sig[i+N_bath] = 1.0*sqrt(w[i]/(2.0*tanh(mu[i])));

    force[0] = F1;         /* assign pointers to force fields */
    force[1] = Fb;
    force[2] = Fb;
    force[3] = F2;
    setwww();




    // ================================================================================================================
    // Processing Path Segments
    // ================================================================================================================

    // Enqueue root path information
    path_info_queue.emplace(PathInfo(-1, 0, 0, 0, 0, random_state));
    // (parent_id -1 for no parent, id, level, clock, random_state)

    // SERIAL IMPLEMENTATION
    // can easily be parallelized:
    // as all paths that are in the queue can be processed in parallel!


    // Loop: all paths breath first
    while (!path_info_queue.empty()) {

        // PARALLEL IMPLEMENTATION
        // [0] -> 0 4 cores: [0,-,-,-]
        // [1,2,3,4] 0 finished -> 1,2,3,4 on 4 cores: [1,2,3,4]
        // [9,10,11,12] 2 finished -> 9 on free core, 4 cores: [1,9,3,4]
        // [10,11,12, 5,6,7,8, 17,18,19,20, 13,14,15,16] 1,4,3 finished -> 10,11,12 on free cores, 4 cores: [10,9,11,12]
        // ...

        // Dequeue path information
        PathInfo path_info = path_info_queue.front();
        path_info_queue.pop();

        // Process path
       /* for (int i = 0; i < Nsample; ++i){
            cout << "Particle:" << i << endl;
            process_path(path_info);
        }*/
        process_path(path_info);

    }



    // ================================================================================================================
    // OUTPUT
    // ================================================================================================================
/*
    long path = n_paths-1; // choose any one between 0 and n_paths-1

    cout << endl;
    cout << "path " << path << endl;

    cout << "valid: " << multi_paths_data[path].valid << endl;
    cout << "parent_id: " << multi_paths_data[path].parent_id << endl;
    for (long i=0; i< multi_paths_data[path].n_data1D; ++i) {
        cout << "multi_paths_data[" << path << "].data1D[" << i << "]: " << multi_paths_data[path].data1D[i] << endl;
    }
    for (long i=0; i< multi_paths_data[path].n_data2D_1; ++i) {
        for (long j=0; j < multi_paths_data[path].n_data2D_2; ++j) {
            cout << "multi_paths_data[" << path << "].data2D[" << i << "][" << j << "]: " << multi_paths_data[path].data2D[i][j] << endl;
        }
    }

    // !!! Multi Paths Data
    // When implemented as C-style dynamic memory (malloc, free) or
    // C++-style dynamic memory (new, delete []),
    // memory needs to be de-allocated explicitly here.
*/
    // =======================================================================================================
    // Memory Deallocation
    // ==========================================================================================================
    delete [] mww;
    delete [] mu;
    delete [] sig;
    delete [] dtdtm;
    delete [] dgam;
    delete [] dhat;
    delete [] R1;
    delete [] v;
    delete [] f;
    delete [] c;
    delete [] m;
    delete [] w;

    return 0;
}
