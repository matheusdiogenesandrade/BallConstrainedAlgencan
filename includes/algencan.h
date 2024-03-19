#ifndef ALGENCAN_H
#define ALGENCAN_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" 
{
#endif
    void c_algencan(void (*myevalf)(int, double*, double*, int*), 
            void (*myevalg)(int, double*, double*, int*), 
            void (*myevalh)(int, double*, int*, int*, double*, int*, int, _Bool*, int*), 
            void (*myevalc)(int, double*, int, double*, int*),
            void (*myevaljac)(int, double*, int, int*, double*, int*, int, _Bool*, int*), 
            void (*myevalhc)(int, double*, int, int*, int*, double*, int*, int, _Bool*, int*), 
            void (*myevalfc)(int, double*, double*, int, double*, int*), 
            void (*myevalgjac)(int, double*, double*, int, int*, int*, double*, int*, int, _Bool*, int*),
            void (*myevalgjacp)(int, double*, double*, int, double*, double*, char, _Bool*, int*), 
            void (*myevalhl)(int, double*, int, double*, double, double*, int*, int*, double*, int*, int, _Bool*, int*),
            void (*myevalhlp)(int, double*, int, double*, double, double*, double*, double*, _Bool*, int*), 
            int jcnnzmax, int hnnzmax, double *epsfeas, double *epsopt, double *efstin, double *eostin, double *efacc, double *eoacc, char *outputfnm, char *specfnm, int nvparam, char **vparam, int n, double *x, double *l, double *u, int m, double *lambda, _Bool *equatn, _Bool *linear, _Bool *coded, _Bool checkder, double *f, double *cnorm, double *snorm, double *nlpsupn, int *inform);

#ifdef __cplusplus
}
#endif

#endif
