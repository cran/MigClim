/*
** migclim.h: Header file for the MigClim methods.
**
** Wim Hordijk   Last modified: 05 October 2011
*/

#ifndef _MIGCLIM_H_
#define _MIGCLIM_H_

/*
** Include files.
*/
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <R.h>


/*
** Defines.
**
** UNIF01:         Draw a uniform random number in [0;1]. Note that 'random'
**                 does not work on Windows %-/  so we use 'rand' instead.
** WEAK_BARRIER:   Weak barrier type.
** STRONG_BARRIER: Strong barrier type.
*/
#define UNIF01         ((double)rand () / RAND_MAX)
#define WEAK_BARRIER   1
#define STRONG_BARRIER 2


/*
** Global variables (we just use many global var's here to avoid passing too
** many arguments all the time).
*/

extern int     nrRows, nrCols, nrEnvChgSteps, nrDispSteps, dispDist, initMatAge,
               fullMatAge, rcThreshold, barrierType, minDist, maxDist;
extern double *dispKernel, *seedProdProb, lddFreq;
extern char    initDistrFile[128], hsMapFile[128], simulName[128],
               barrierFile[128];
extern bool    useBarrier, fullOutput;


/*
** Function prototypes.
*/
void mcMigrate           (char **paramFile, int *nrFiles);
bool mcSrcCell           (int i, int j, int **curState, int **pxlAge,
			  int loopID, int habSuit, int **barriers);
int  mcUnivDispCnt       (int **habSuit, int **barriers);
void updateNoDispMat     (int **hsMat, int **noDispMat, int *noDispCount);
void mcFilterByBarrier   (int **curState, int **barriers);
bool mcIntersectsBarrier (int snkX, int snkY, int srcX, int srcY,
			  int **barriers);
int  mcInit              (char *paramFile);
int  mcReadMatrix        (char *fname, int **mat);
int  mcWriteMatrix       (char *fname, int **mat);


#endif  /* _MIGCLIM_H_ */

/*
** EoF: migclim.h
*/
