#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <complex>


#include "Global.h"
#include "RandomNumberGenerator.h"
#include "MultiPathProcessing.h"
#include "variable.h"
#include "variable-trajectory.h"

#include "dmatrix.c"
#include "imatrix.c"

#include "density.cpp"
#include "functions.cpp"
#include "transition.cpp"


// Global Data
extern PathInfoQueue_t  path_info_queue;
extern PathDataVector_t multi_paths_data;

using namespace std;

#include <gsl/gsl_rng.h>


char  datafilename[80];
FILE   *stream;

//double  *abszsum1, *argzsum1; //0 - summand trajectory dependent ; 1 - sum
//double  *habszsum1, *hargzsum1;

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

    // !!! INPUT

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
    Nsample = 10000;
    w_max = 3;
    eta = 0.13;
    beta = 25;
    delta = 0.8;
    ppower = 100000;


    // Random Number Generator
    // fixed seed for reproducibility, otherwise use RandomState()
    cout << "Root Seed: " << init_seed << endl;
    RandomState random_state = RandomState(init_seed); // root path

    // !!! Multi Paths Data
    long n_data1D   = 4; // dimension parameters // inital values needed by children
    long n_data2D_1 = N_slice; // sum1 length
    long n_data2D_2 = 4; // abszum, argzum, habszum, hargzsum
    long n_paths = pow(N_PATHS, (N_LEVELS+1.0)) - 1;

    multi_paths_data.resize(n_paths, PathData(n_data1D, n_data2D_1, n_data2D_2));

    // ======================= Memory Allocation ================================




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
        process_path(path_info);
    }

    // !!! OUTPUT

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


    // ======================== Memory Deallocation ================================



    return 0;
}
