#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <vector>
#include <iostream>

#include "algencan.h"

using namespace std;

// data
int N = 12;
   
double Distance[] = {
        0, 1, 2, 2, 3, 4, 4, 5, 3, 5, 6, 7,
        1, 0, 1, 1, 2, 3, 3, 4, 2, 4, 5, 6,
        2, 1, 0, 2, 1, 2, 2, 3, 1, 3, 4, 5,
        2, 1, 2, 0, 1, 2, 2, 3, 3, 3, 4, 5,
        3, 2, 1, 1, 0, 1, 1, 2, 2, 2, 3, 4,
        4, 3, 2, 2, 1, 0, 2, 3, 3, 1, 2, 3,
        4, 3, 2, 2, 1, 2, 0, 1, 3, 1, 2, 3,
        5, 4, 3, 3, 2, 3, 1, 0, 4, 2, 1, 2,
        3, 2, 1, 3, 2, 3, 3, 4, 0, 4, 5, 6,
        5, 4, 3, 3, 2, 1, 1, 2, 4, 0, 1, 2,
        6, 5, 4, 4, 3, 2, 2, 1, 5, 1, 0, 1,
        7, 6, 5, 5, 4, 3, 3, 2, 6, 2, 1, 0
};
   
double Flow[] = {
    0, 3, 4, 6, 8, 5, 6, 6, 5, 1, 4, 6,
    3, 0, 6, 3, 7, 9, 9, 2, 2, 7, 4, 7,
    4, 6, 0, 2, 6, 4, 4, 4, 2, 6, 3, 6,
    6, 3, 2, 0, 5, 5, 3, 3, 9, 4, 3, 6,
    8, 7, 6, 5, 0, 4, 3, 4, 5, 7, 6, 7,
    5, 9, 4, 5, 4, 0, 8, 5, 5, 5, 7, 5,
    6, 9, 4, 3, 3, 8, 0, 6, 8, 4, 6, 7,
    6, 2, 4, 3, 4, 5, 6, 0, 1, 5, 5, 3,
    5, 2, 2, 9, 5, 5, 8, 1, 0, 4, 5, 2,
    1, 7, 6, 4, 7, 5, 4, 5, 4, 0, 7, 7,
    4, 4, 3, 3, 6, 7, 6, 5, 5, 7, 0, 9,
    6, 7, 6, 6, 7, 5, 7, 3, 2, 7, 9, 0
};


// helpers
int getIndex(int s, int i) {
    return s * N + i; 
}

double getDistance(int s, int i) {
    return Distance[getIndex(s, i)];
}

double getFlow(int s, int i) {
    return Flow[getIndex(s, i)];
}

double& getElement(double * array, int s, int i) {
    return array[getIndex(s, i)];
}

int getVar(double * x, int s, int i) {
    return (x[getIndex(s, i)] + 1.)/2.;
}

bool enoughMemory(int * nnz, int lim, bool * lmem) {
    if( *nnz > lim ) {
        *lmem = 1;
        return false;
    }
    return true;
}

/* ******************************************************************
 *  ****************************************************************** */

void myevalf(int n, double *x, double *f, int *flag) {

    *flag = 0;

    *f = 0;
    for (int s = 0; s < N; ++s)
        for (int t = 0; t < N; ++t)
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    *f += Distance[getIndex(i, j)] * Flow[getIndex(s, t)] * getVar(x, s, i) * getVar(x, t, j);
}

/* ******************************************************************
 *  ****************************************************************** */

void myevalg(int n, double *x, double *g, int *flag) {

    *flag = 0;

    for (int s = 0; s < N; ++s) 
        for (int i = 0; i < N; ++i) {
            double& entry = getElement(g, s, i);

            for (int t = 0; t < N; ++t)
                for (int j = 0; j < N; ++j) 
                    entry += (2. * getVar(x, t, j)) * getDistance(i, j) * getFlow(s, t);

        }
}

/* ******************************************************************
 *  ****************************************************************** */

void myevalh(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz,
        int lim, _Bool *lmem, int *flag) {

    *flag = 0;
    *lmem = 0;
    *hnnz = (N * N * N * N) / 2.0;

    // edge case
    if (!enoughMemory(hnnz, lim, lmem)) return;

    // fill 
    int hessian_idx = 0;
    for (int s = 0; s < N; ++s) 
        for (int t = 0; t < N; ++t)
            for (int i = 0; i < N; ++i) 
                for (int j = 0; j < N; ++j) {
                    int row = getIndex(s, i);
                    int col = getIndex(t, j);

                    if (row >= col) {
                        hrow[hessian_idx] = row;
                        hcol[hessian_idx] = col;

                        hval[hessian_idx] = getDistance(i, j) * getFlow(s, t);

                        hessian_idx++;
                    }
                }
}

/* ******************************************************************
 *  ****************************************************************** */

void myevalc(int n, double *x, int ind, double *c, int *flag) {

    *flag = 0;

    // fill 
    if (ind < N) {
        /* FACILITIES CONSTRAINTS */
    
        // get facility
        const int s = ind;

        // set value
        *c = 1;
        for (int i = 0; i < N; ++i)
            *c -= getVar(x, s, i);
        // 1 - \sum_{i = 1,...,n} (x_{si} + 1)/2 = 0

    } else if (ind < 2 * N) {
        /* LOCATIONS CONSTRAINTS */

        // get location
        const int i = ind % N;

        // set value
        *c = 1;
        for (size_t s = 0; s < N; ++s)
            *c -= getVar(x, s, i);
        // 1 - \sum_{s = 1,...,n} (x_{si} + 1)/2 = 0

    } else if (ind < 2 * N + 1) {
        /* Ball constraint */
    
        // set value
        *c = N * N;
        double sum = 0;
        for (int s = 0; s < N; ++s)
            for (int i = 0; i < N; ++i)
                *c -= getElement(x, s, i) * getElement(x, s, i);
        // n * n - \sum_{s, i = 1,...,n} x_{si} = 0

    } else {
        // invalid index
        *flag = -1;
    }
}

/* ******************************************************************
 *  ****************************************************************** */

void myevaljac(int n, double *x, int ind, int *jcvar, double *jcval,
        int *jcnnz, int lim, _Bool *lmem, int *flag) {

    *flag = 0;
    *lmem = 0;

    // fill 
    if (ind < N) {
        /* FACILITIES CONSTRAINTS */

        // get facility
        const int s = ind;

        // edge case
        *jcnnz = N;
        if (!enoughMemory(jcnnz, lim, lmem)) return;

        // set value
        for (int i = 0; i < N; ++i) {
            jcvar[i] = getIndex(s, i);
            jcval[i] = - 1. / 2.;
        }

    } else if (ind < 2 * N) {
        /* LOCATIONS CONSTRAINTS */
    
        // get facility
        const int i = ind % N;

        // edge case
        *jcnnz = N;
        if (!enoughMemory(jcnnz, lim, lmem)) return;

        // set value
        for (int s = 0; s < N; ++s) {
            jcvar[s] = getIndex(s, i);
            jcval[s] = - 1. / 2.;
        }

    } else if (ind < 2 * N + 1) {
        /* Ball constraint */
    
        // edge case
        *jcnnz = N * N;
        if (!enoughMemory(jcnnz, lim, lmem)) return;

        // set value
        size_t idx = 0;
        for (size_t s = 0; s < N; ++s)
            for (size_t i = 0; i < N; ++i) {
                jcvar[idx] = getIndex(s, i);
                jcval[idx] = - 2. * getElement(x, s, i);
                ++idx;
            }

    } else {
        // invalid index
        *jcnnz = 0;
//        *flag = -1;
    }
}

/* ******************************************************************
 *  ****************************************************************** */

void myevalhc(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval,
        int *hcnnz, int lim, _Bool *lmem, int *flag) {

    *flag = 0;
    *lmem = 0;

    // edge case 1
    if (ind != 2 * N) {
        *hcnnz = 0;
        return;
    }

    // edge case 2
    *hcnnz = N * N;
    if (!enoughMemory(hcnnz, lim, lmem)) return;

    // fill 
    int hessian_idx = 0;
    for (int s = 0; s < N; ++s) 
        for (int i = 0; i < N; ++i) {
            int var_idx = getIndex(s, i);

            hccol[hessian_idx] = var_idx;
            hcrow[hessian_idx] = var_idx;
            hcval[hessian_idx] = - 1.;
            hessian_idx++;
        }
}

/* *****************************************************************
 *  ***************************************************************** */

void myevalfc(int n, double *x, double *f, int m, double *c, int *flag) { *flag = -1; }

/* *****************************************************************
 *  ***************************************************************** */

void myevalgjac(int n, double *x, double *g, int m, int *jcfun, int *jcvar,
        double *jcval, int *jcnnz, int lim, _Bool *lmem, int *flag) { *flag = -1; }

/* *****************************************************************
 *  ***************************************************************** */

void myevalgjacp(int n, double *x, double *g, int m, double *p, double *q,
        char work, _Bool *gotj, int *flag) { *flag = -1; }

/* *****************************************************************
 *  ***************************************************************** */

void myevalhl(int n, double *x, int m, double *lambda, double scalef,
        double *scalec, int *hlrow, int *hlcol, double *hlval,
        int *hlnnz, int lim, _Bool *lmem, int *flag) { *flag = -1; }

/* *****************************************************************
 *  ***************************************************************** */

void myevalhlp(int n, double *x, int m, double *lambda, double scalef,
        double *scalec, double *p, double *hp, _Bool *goth, 
        int *flag) { *flag = -1; }

/* ******************************************************************
 *  ****************************************************************** */

int main() {
    //
    const int m = 2 * N + 1;
    const int n = N * N;


    /* Memory allocation */
    vector<double> x (n, -1.);
    vector<double> l (n, -1.);
    vector<double> u (n, 1.);
    vector<double> lambda (m, 0.);
    bool * equatn = (bool *) calloc(m, sizeof(bool));
    bool * linear = (bool *) calloc(m, sizeof(bool));

    /* Initial point */
    vector<int> seq = {3,10,11,2,12,5,6,7,8,1,4,9};
    for (int i = 0; i < seq.size(); ++i)
        getElement(x.data(), seq[i] - 1, i) = 1.;

    /* Initial lambda */
    vector<double> lagrangian_coeffs = {
        504.0,
        292.0,
        192.0,
        232.0,
        300.0,
        364.0,
        184.0,
        420.0,
        564.0,
        572.0,
        484.0,
        660.0,
        0.0,
        204.0,
        220.0,
        -112.0,
        220.0,
        264.0,
        296.0,
        8.0,
        32.0,
        168.0,
        100.0,
        256.0,
        -21.0
    };
    for (int i = 0; i < m; ++i)
        lambda[i] = lagrangian_coeffs[i];

    //
    for (int k = 0; k < m; ++k)
        equatn[k] = true;

    //
    for (int k = 0; k < m - 1; ++k)
        linear[k] = true;

    //
    bool  coded[11];
    coded[0]  = 1; /* fsub     */
    coded[1]  = 1; /* gsub     */
    coded[2]  = 1; /* hsub     */
    coded[3]  = 1; /* csub     */
    coded[4]  = 1; /* jacsub   */
    coded[5]  = 1; /* hcsub    */
    coded[6]  = 0; /* fcsub    */
    coded[7]  = 0; /* gjacsub  */
    coded[8]  = 0; /* gjacpsub */
    coded[9]  = 0; /* hlsub    */
    coded[10] = 0; /* hlpsub   */

    //
    int jcnnzmax = 4 * N * N;
    int hnnzmax  = N * N + (N * N * N * N + N * N) / 2.;

    //
    bool checkder = false;

    /* Parameters setting */
    double epsfeas = 1.0e-08;
    double epsopt  = 1.0e-08;
    double efstin  = sqrt( epsfeas );
    double eostin  = pow( epsopt, 1.5 );
    double efacc   = sqrt( epsfeas );
    double eoacc   = sqrt( epsopt );

    //
    char *outputfnm = "algencan.out";
    char *specfnm = "";

    int nvparam = 2;
//    int nvparam = 17;

    //
    char **vparam = ( char ** ) malloc( nvparam * sizeof( char * ) );
    vparam[0]  = "ITERATIONS-OUTPUT-DETAIL 10";
    vparam[1]  = "OBJECTIVE-AND-CONSTRAINTS-SCALING-AVOIDED";

//    vparam[1]  = "INNER-ITERATIONS-LIMIT 1000";
//    vparam[2]  = "OUTER-ITERATIONS-LIMIT 1000";
//    vparam[3]  = "NUMBER-OF-ARRAYS-COMPONENTS-IN-OUTPUT 6";
//    vparam[4]  = "ACCELERATION-PROCESS-ITERATIONS-LIMIT 1000";
//    vparam[5]  = "PENALTY-PARAMETER-INITIAL-VALUE 1e-5";
//    vparam[6]  = "LARGEST-PENALTY-PARAMETER-ALLOWED 1e22";
//    vparam[7]  = "SOLUTION-FILENAME sol.sol";
//    vparam[8]  = "OBJECTIVE-AND-CONSTRAINTS-SCALING-AVOIDED";
//    vparam[9]  = "FIXED-VARIABLES-REMOVAL-AVOIDED";
//    vparam[10] = "LINEAR-SYSTEMS-SOLVER-IN-ACCELERATION-PROCESS MA97";
//    vparam[11] = "TRUST-REGIONS-INNER-SOLVER MA57";
//    vparam[12] = "LINEAR-SYSTEMS-SOLVER-IN-TRUST-REGIONS MA97";
//    vparam[13] = "NEWTON-LINE-SEARCH-INNER-SOLVER MA97";
//    vparam[14] = "LINEAR-SYSTEMS-SOLVER-IN-NEWTON-LINE-SEARCH MA97";
//    vparam[15] = "TRUNCATED-NEWTON-LINE-SEARCH-INNER-SOLVER TRUEHP";
//    vparam[16] = "MATRIX-VECTOR-PRODUCT-IN-TRUNCATED-NEWTON-LS TRUEHP";

//    vparam[13]  = "SKIP-ACCELERATION-PROCESS";
//    vparam[15]  = "ADD-SLACKS";
//    vparam[16]  = "IGNORE-OBJECTIVE-FUNCTION,
    //
    double cnorm,f,nlpsupn,snorm;
    int inform;

    //
    c_algencan(
            &myevalf,
            &myevalg,
            &myevalh,
            &myevalc,
            &myevaljac,
            &myevalhc,
            &myevalfc,
            &myevalgjac,
            &myevalgjacp,
            &myevalhl,
            &myevalhlp,
            jcnnzmax,
            hnnzmax,
            &epsfeas,
            &epsopt,
            &efstin,
            &eostin,
            &efacc,
            &eoacc,
            outputfnm,
            specfnm,
            nvparam,
            vparam,
            n,
            x.data(),
            l.data(),
            u.data(),
            m,
            lambda.data(),
            equatn,
            linear,
            coded,
            checkder,
            &f,
            &cnorm,
            &snorm,
            &nlpsupn,
            &inform
                );

    /* Memory deallocation */
    free(equatn);
    free(linear);
}

