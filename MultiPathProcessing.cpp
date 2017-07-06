#include "MultiPathProcessing.h"
#include <complex>

// Global Data
PathInfoQueue_t  path_info_queue;
PathDataVector_t multi_paths_data;

// =============================================================================
// Path Processing
// =============================================================================

void process_path(PathInfo& path_info) {

    cout << endl;
    cout << "Process Path " << path_info.id << " (level " << path_info.level << " )" << endl;
    cout << "No. of jumps in path: " << path_info.Njump << endl;

    // Read ancestors PathData
    long ancestor_id = path_info.parent_id;
    while (ancestor_id >= 0) { // root path id is 0, parent of root path id is -1
        cout << "ancestor: " << ancestor_id << endl;
        // How to access PathData
        // multi_paths_data[ancestor_id].valid ... always valid
        // multi_paths_data[ancestor_id].parent_id
        // multi_paths_data[ancestor_id].data1D[i]
        // multi_paths_data[ancestor_id].data2D[i][j]
        ancestor_id = multi_paths_data[ancestor_id].parent_id; // next ancestor
    }

    // !!! PROCESSING
    complex <double> z;
    /* Jump condition to be placed here;
     *  //calculating matrix elements

        //avoid index problem if X < # then = 0, otherwise Y
        ap0 = fabs(p0 = ((www[1][SS0][0]() < -7775.0) ? 0.0 : www[0][SS0][0]()));
        ap1 = fabs(p1 = ((www[1][SS0][1]() < -7775.0) ? 0.0 : www[0][SS0][1]()));
        ap2 = fabs(p2 = ((www[1][SS0][2]() < -7775.0) ? 0.0 : www[0][SS0][2]()));
        ap3 = fabs(p3 = ((www[1][SS0][3]() < -7775.0) ? 0.0 : www[0][SS0][3]()));
        dn2 = ap0 + ap1 + ap2 + ap3;

        xx = dn2 * (gsl_rng_uniform(rr));   // choosing matrix elements
        // printf(" oldz  real %lf imag %lf\t", real(z), imag(z));

        //similar to perturbation expansion
        //keeping track of jumps
        //corralted sapling needed here
        //for certain x will get unlikely jumps
        oldz = z;
        SS2 = SS1;
        if (xx < ap0) {
            SS1 = 0;
            z *= p0 * dn2 / ap0;
        } else if (xx < ap0 + ap1) {
            SS1 = 1;
            z *= p1 * dn2 / ap1;
        } else if (xx < ap0 + ap1 + ap2) {
            SS1 = 2;
            z *= p2 * dn2 / ap2;
        } else {
            SS1 = 3;
            z *= p3 * dn2 / ap3;
        }
        if (SS0 != SS1)
            NNjmp++;
        if (NNjmp > Ncut)
            return 0;

      My major issue is I'm not sure how this is implemented. Discuss with Athina, once other modules are added then it should be a simple
      call job.

      There is a need for the child production. Here's the idea;

      if (SS0 != SS1)
            NNjmp++ and produce 4 children;
      if (NNjmp > Ncut)
            return 0;

        */

    long clock = path_info.clock;
    while (clock < N_CLOCKS) {

        double r = path_info.random_state.uniform_real(0.0, 1.0);
        cout << "clock: " << clock << " random number: " << r << endl;
        clock++;
        if (r > 0.95) break;
    }

    // Write PathData
    multi_paths_data[path_info.id].valid = true;
    multi_paths_data[path_info.id].parent_id = path_info.parent_id;
    multi_paths_data[path_info.id].probability = z;
    for (long i=0; i< multi_paths_data[path_info.id].n_data1D; ++i) {
        multi_paths_data[path_info.id].data1D[i] = path_info.id;
    }
    for (long i=0; i< multi_paths_data[path_info.id].n_data2D_1; ++i) {
        for (long j=0; j < multi_paths_data[path_info.id].n_data2D_2; ++j) {
            multi_paths_data[path_info.id].data2D[i][j] = path_info.id;
        }
    }

    if ((path_info.Njump < N_JUMPS) &&  (clock < N_CLOCKS)) {

        // Calculate the following paths ids using a formula (vs. using a shared counter)
        // to avoid synchronization between parallel executions of process_path().
        // lowest path id in current level
        long id_min_level = 0;
        for (long l=0; l<path_info.level; ++l) id_min_level += pow(N_PATHS, l);
        // lowest path id in next level
        long id_min_next_level = id_min_level + pow(N_PATHS, path_info.level);
        // first path id of following paths in next level
        long path_id = id_min_next_level + (path_info.id - id_min_level)*N_PATHS;

        // Enqueue path following paths information
        for (long p=0; p<(N_PATHS-1); ++p) {
            // Random Number Generator
            // new random states based on different seeds create with random state from current path
            unsigned long seed = path_info.random_state.uniform_int(0, path_info.random_state.MAX_INT);
            cout << "Child Seed " << p << ": " << seed << endl;
            RandomState temp_random_state = RandomState(seed);
            path_info_queue.emplace(PathInfo(path_info.id, path_id + p, path_info.level + 1, path_info.Njump + 1,  clock, temp_random_state));
            // (parent_id, id, level, clock, random_state)
        }
        // Pass on random state from current path to one of the following paths
        // after the random state has been used to generate new seeds for the following paths
        // to make sure that different seeds are generated in the following path.
        path_info_queue.emplace(PathInfo(path_info.id, path_id + (N_PATHS-1), path_info.level+1, path_info.Njump, clock, path_info.random_state));
        // (parent_id, id, level, clock, random_state)
    }

    return;
}