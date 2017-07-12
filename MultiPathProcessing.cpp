#include "MultiPathProcessing.h"

using namespace std;

// Global Data
PathInfoQueue_t  path_info_queue;
PathDataVector_t multi_paths_data;

// =============================================================================
// Variables
// =============================================================================

// Local =================
int SS0;
int SS1;
int SS3;

double phase0 = 0.0;
double alpha = 0.0;
double p0;
double p1;
double p2;
double p3;
double ap0;
double ap1;
double ap2;
double ap3;
double dn2;
double xx;
double signPdotdhat;
complex<double> oldz;

double *abszsum1;
double *argzsum1;
double *habszsum1;
double *hargzsum1;
double *Pperp;

double (*phi)(double *, double *, int);

// Extern ================


extern int N_bath;
extern int N_slice;
extern double sina;
extern double cosa;

extern double Dt;
extern double de;
extern double Pdotdhat;
extern double abs_d;

extern double *R1;
extern double *v;
extern double *dhat;

extern const gsl_rng_type * TT;
extern gsl_rng * rr;

extern complex<double> z;
extern complex<double> I;

extern double (*dens_init[4])(double*, double*, int );
extern double (*obs[4])(double*, double*, int );
extern double (*obs1[4])(double*, double*, int );

extern double (* www[6][4][4])();



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


    abszsum1 = new double[N_slice];
    argzsum1 = new double[N_slice];
    habszsum1 = new double[N_slice];
    hargzsum1 = new double[N_slice];
    Pperp = new double[N_bath];

    for (int i = 0; i < N_slice; ++i){
        abszsum1[i] = 0.0;
        argzsum1[i]  = 0.0;
        habszsum1[i] = 0.0;
        hargzsum1[i] = 0.0;
    }

    // ==============================================================================================================
    // Initial Distribution
    // ==============================================================================================================
    if (path_info.id == 0){
        gauss_init_W(R1, v);

        double yy = 4.0 * (gsl_rng_uniform(rr));
        if (yy < 1.0)
            SS3 = (SS0 = 0);
        else if (yy < 2.0) {
            SS0 = 1;
            SS3 = 2;
        } else if (yy < 3.0) {
            SS0 = 2;
            SS3 = 1;
        } else
            SS3 = (SS0 = 3);

        z = 4.0;
        SS1 = SS0;
        multi_paths_data[path_info.id].surface = SS0;
        cout << "Initial Surface: " << SS1 <<endl;

    }


    // =======================================================================================
    // JUMP CONDITION
    // =======================================================================================
    long clock = path_info.clock;
    int N_CLOCKS = N_slice;
    do{
        SS0 = multi_paths_data[path_info.id].surface;
        cout << "New original surface: " << SS0 << endl;
        phase0 = U(R1, v, SS0, Dt * 0.5);
        z *= exp(I * phase0);

        dd(R1); //0 corresponds to the # of ensembles (single path) - needs to be defined for every path
        de = dE(R1);

        for (int i = 0, Pdotdhat = 0; i < N_bath; ++i) {
            Pdotdhat += v[i] * dhat[i];
        }
        alpha = Pdotdhat * abs_d * Dt;


        signPdotdhat = (Pdotdhat < 0 ? -1 : 1); // -1 if Pdotdhat neg, 1 is pos
        Pdotdhat = fabs(Pdotdhat);
        for (int i = 0; i < N_bath; ++i)
            Pperp[i] = v[i] - signPdotdhat * Pdotdhat * dhat[i];

        alpha *= 2.0;
        sina = sin(alpha);
        cosa = cos(alpha);

        //calculating matrix elements

        //avoid index problem if X < # then = 0, otherwise Y
        ap0 = fabs(p0 = ((www[1][SS0][0]() < -7775.0) ? 0.0 : www[0][SS0][0]()));
        ap1 = fabs(p1 = ((www[1][SS0][1]() < -7775.0) ? 0.0 : www[0][SS0][1]()));
        ap2 = fabs(p2 = ((www[1][SS0][2]() < -7775.0) ? 0.0 : www[0][SS0][2]()));
        ap3 = fabs(p3 = ((www[1][SS0][3]() < -7775.0) ? 0.0 : www[0][SS0][3]()));
        dn2 = ap0 + ap1 + ap2 + ap3;
    //    cout << "Prob:" << "ap0" << ap0 <<" ap1"<< ap1 <<" ap2" << ap2 <<" ap3"<< ap3 << endl;
        xx = dn2 * (gsl_rng_uniform(rr));   // choosing matrix elements

        //similar to perturbation expansion
        //keeping track of jumps
        oldz = z;
        if (xx < ap0) {
            SS1 = 0;
            z *= p0 * dn2 / ap0;
            cout << "New Path: "<< SS1 << endl;
        } else if (xx < ap0 + ap1) {
            SS1 = 1;
            z *= p1 * dn2 / ap1;
            cout << "New Path: "<< SS1 << endl;
        } else if (xx < ap0 + ap1 + ap2) {
            SS1 = 2;
            z *= p2 * dn2 / ap2;
            cout << "New Path: "<< SS1 << endl;
        } else {
            SS1 = 3;
            z *= p3 * dn2 / ap3;
            cout << "New Path: "<< SS1 << endl;
        }
        /*if (www[1][SS0][SS1]() != 9999.0)
            for (int i = 0; i < N_bath; ++i)
                v[i] = Pperp[i] + signPdotdhat * www[1][SS0][SS1]() * dhat[i];
         */

        // ================================================================================
        // Cummulative Sum - done for each timeslice
        // ================================================================================

        double abszsum0;
        double argzsum0;
        double habszsum0;
        double hargzsum0;

        phase0 = U(R1, v, SS1, Dt*0.5);
        z *= exp(I * phase0);
        phi = obs[SS1];

        abszsum0 = real(z * phi(R1, v, 0) * dens_init[SS3](R1, v, 1));
        argzsum0 = imag(z * phi(R1, v, 0) * dens_init[SS3](R1, v, 1));
        abszsum1[clock] += abszsum0;
        argzsum1[clock] += argzsum0;

        phi = obs1[SS1];
        habszsum0 = real(z * phi(R1, v, 0) * dens_init[SS3](R1, v, 1));
        hargzsum0 = imag(z * phi(R1, v, 0) * dens_init[SS3](R1, v, 1));
        habszsum1[clock] += habszsum0;
        hargzsum1[clock] += hargzsum0;


        multi_paths_data[path_info.id].surface = SS1;
        clock++;
    } while (clock < N_slice || SS1 != SS0);


/*
    while (clock < N_CLOCKS) {

        double r = path_info.random_state.uniform_real(0.0, 1.0);
        cout << "clock: " << clock << " random number: " << r << endl;

        clock++;
        if (r > 0.95) break;
    }
*/

    // Write PathData
    multi_paths_data[path_info.id].valid = true;
    multi_paths_data[path_info.id].parent_id = path_info.parent_id;
    multi_paths_data[path_info.id].probability = z;
    multi_paths_data[path_info.id].phase = phase0;
    multi_paths_data[path_info.id].surface = SS1;
    for (long i = 0; i < multi_paths_data[path_info.id].n_data1D; ++i) { //initial density values for children
        multi_paths_data[path_info.id].data1D[i] = path_info.id;
    }
    for (long i = 0; i < multi_paths_data[path_info.id].n_data2D_2; ++i) {
            multi_paths_data[path_info.id].data2D[0][i] = abszsum1[i];
            multi_paths_data[path_info.id].data2D[1][i] = argzsum1[i];
            multi_paths_data[path_info.id].data2D[2][i] = habszsum1[i];
            multi_paths_data[path_info.id].data2D[3][i] = hargzsum1[i];
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

    delete [] abszsum1; delete [] argzsum1; delete [] habszsum1; delete [] hargzsum1; delete [] Pperp;

    return;
}