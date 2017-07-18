#ifndef PATHPROCESSING_H
#define PATHPROCESSING_H

#include <iostream>
#include <vector>
#include <queue>
#include <complex>
#include <gsl/gsl_rng.h>

#include "Global.h"
#include "RandomNumberGenerator.h"

using namespace std;

//extern double  *abszsum1, *argzsum1; //0 - summand trajectory dependent ; 1 - sum
//extern double  *habszsum1, *hargzsum1;

extern int N_slice;

// =============================================================================
// Path Processing Queue (FIFO)
// =============================================================================

// A queue (FIFO) provides breath first processing of the path structure (tree),
// vs. a stack (LIFO) that provides depth first processing.
// Breath first processing is better for parallel processing,
// as it adds more paths to the processing queue quicker.

struct PathInfo {
    PathInfo(long parent_id, long id, long surface, complex<double> z, long level, int Njump, long clock/*, RandomState random_state*/)
            : parent_id(parent_id), id(id), surface(surface), probability(z), level(level), Njump(Njump), clock(clock)/*, random_state(random_state)*/ { };
    long parent_id;
    long id;
    long surface;
    complex<double> probability;
    long level;
    int Njump;
    long clock;
    RandomState random_state;
};

typedef queue<PathInfo>  PathInfoQueue_t;

// =============================================================================
// !!! Multi Paths Data
// =============================================================================

// Shared memory for process_path()for writing and reading final and shared results.

// No data protection is needed when parallelizing process_path(),
// because only one path is writing the data, and only following paths,
// that do not exist yet when the data is written, will only read the data.

// Can be implemented as C-style dynamic memory (malloc, free),
// C++-style dynamic memory (new, delete []), unique pointers (C++11), or
// STL containers, like vectors, arrays (C++11), lists.

// This example of a vector implementation assumes
// that the PathData for each path has the same dimension.

struct PathData {

    PathData(long n_data2D_1, long n_data2D_2) // dimension parameters of the PathData
            :n_data2D_1(n_data2D_1), n_data2D_2(n_data2D_2),
              valid(false), parent_id(-1), probability(4, (1.0, 0.0)), surface(4, 0),
                abszsum1(n_data2D_1,vector<double>(n_data2D_2, 0.0)),
                argzsum1(n_data2D_1,vector<double>(n_data2D_2, 0.0)),
                habszsum1(n_data2D_1,vector<double>(n_data2D_2, 0.0)),
                hargzsum1(n_data2D_1,vector<double>(n_data2D_2, 0.0)) { }

    PathData(long n_data1D, long n_data2D_1, long n_data2D_2) // dimension parameters of the PathData
            : n_data1D(n_data1D), n_data2D_1(n_data2D_1), n_data2D_2(n_data2D_2),
              valid(false), parent_id(-1), probability(1.0, 0.0), data1D(n_data1D, 0.0), data2D(n_data2D_1,vector<double>(n_data2D_2, 0.0)) { }
  
    // memory is allocated using constructor initialization
    // explicit initialization of vectors with 0.0
    bool                    valid; // identify processed paths
    // given that N_CLOCKS is a stop criteria, not all levels might be reached
    long                    parent_id;
    vector<complex<double>> probability;
    vector<double>          surface;

    vector<vector<double>>  abszsum1;
    vector<vector<double>>  argzsum1;
    vector<vector<double>>  habszsum1;
    vector<vector<double>>  hargzsum1;
    // dimension parameters of the PathData
    long                    n_data2D_1;
    long                    n_data2D_2;
};

typedef vector<PathData> PathDataVector_t;

// =============================================================================
// Path Processing
// =============================================================================

void process_path(PathInfo& path_info);

#endif
