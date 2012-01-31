/*
** migrate.c: Implementation of the main MigClim function.
**
** Wim Hordijk   Last modified: 24 January 2012
**
** This C code is based on the original Visual Basic code of Robin Engler.
*/

#include "migclim.h"


/*
** Global variables, so we won't have to pass too many arguments all the time.
*/
int     nrRows, nrCols, nrEnvChgSteps, nrDispSteps, dispDist, initMatAge,
        fullMatAge, rcThreshold, barrierType, minDist, maxDist;
double *dispKernel, *seedProdProb, lddFreq;
char    initDistrFile[128], hsMapFile[128], simulName[128], barrierFile[128];
bool    useBarrier, fullOutput;
typedef struct _pixel
{
  int row, col;
} pixel;
pixel rndPixel;


/*
** Function prototypes.
*/
void mcRandomPixel   (pixel *pix);
bool mcSinkCellCheck (pixel pix, int **curState, int **habSuit, int **bars);


/*
** mcMigrate: The core of the MigClim method. Perform the main migration steps.
**            Parameter values are read from a file.
**
** Parameters:
**   - paramFile: The name of the parameter file.
**   - nrFiles:   A pointer to an integer to contain the number of output
**                files created. A value of -1 is returned if an error occurred.
*/

void mcMigrate (char **paramFile, int *nrFiles)
{
  int     i, j, envChgStep, dispStep, loopID, simulTime;
  bool    advOutput, habIsSuitable, cellInDispDist, tempResilience;
  char    fileName[128];
  FILE   *fp=NULL, *fp2=NULL;
  double  lddSeedProb;
  time_t  startTime;
  /*
  ** These variables are not (yet) used.
  **
  ** int vegResTime, seedBankResTime, *seedBank, tSinceDecol, tWhenCol;
  */
  
  /*
  ** Pixel counter variables. These variables allow us to record interesting
  ** values to be printed into the output file:
  **
  ** ==> Many of these aren't used yet. <==
  **
  **   - nrColonized:           The number of pixels in "Colonized" state at
  **                            any given time.
  **   - nrAbsent:              The number of pixels in "Absent" state at any
  **                            given time.
  **   - nrDecolonized:         The number of pixels in "Decolonized" state at
  **                            any given time.
  **   - nrTmpResilient:        The number of pixels that in "Temporary
  **                            Resilience" state at any given time.
  **   - nrVegResilient:        The number of pixels that are in "Vegetative
  **                            Resilience" state at any given time.
  **   - nrSeedBank:            The number of pixels that are in Seed Bank
  **                            Resilience state at any given time.
  **   - nrStepColonized:       The number of pixels which changed to
  **                            "Colonized" state during a given step.
  **   - nrStepDecolonized:     The number of pixels which changed to
  **                            "Decolonized" state during a given step.
  **   - nrStepTmpResilient:    The number of pixels that which changed to
  **                            "Temporary Resilience" state during a given
  **                            step.
  **   - nrStepVegResilient:    The number of pixels that are which changed to
  **                            "Vegetative Resilience" state during a given
  **                            step.
  **   - nrStepSeedBank:        The number of pixels that are which changed to
  **                            "Seed Bank Resilience" state during a given
  **                            step.
  **   - nrStepLostUnsuit:      The number of pixels that were lost due to
  **                            unsuitable habitat conditons during a given
  **                            dispersal step.
  **   - nrStepLostExtreme:     The number of pixels that were lost due to
  **                            "Extreme Event" during a given dispersal step.
  **   - nrStepLDDSuccess:      The number of LDD events that were successful
  **                            (can be any of a "normal" LDD or a River or
  **                            Road LDD).
  **   - nrStepVegResRecover:   The number of "Vegetative Resilience" recovery
  **                            events that occured within a given dispersal
  **                            step.
  **   - nrStepSeedBankRecover: The number of "Seed Bank Resilience" recovery
  **                            events that occured within a given dispersal
  **                            step.
  **   - nrTotColonized:        \
  **   - nrTotLostUnsuit:        |
  **   - nrTotLostExtreme:       |   Same as above but over the total
  **   - nrTotLostNatural:        >  simulation instead of one step.
  **   - nrTotLDDSuccess:        |
  **   - nrTotVegResRecover:     |
  **   - nrTotSeedBankRecover:  /
  **   - nrInitial:             The number of initial pixels that are occupied
  **                            by the species.
  **   - nrNoDispersal:         The number of pixels that would be colonized at
  **                            the end of the simulation under the
  **                            "no-dispersal" hypothesis.
  **   - nrUnivDispersal:       The number of pixels that would be colonized
  **                            at the end of the simulation under the
  **                            "unlimited-dispersal" hypothesis.
  */
  int nrColonized, nrAbsent, nrStepColonized,
      nrStepDecolonized, nrStepLDDSuccess, 
      nrTotColonized, nrTotLDDSuccess, nrInitial, nrNoDispersal,
      nrUnivDispersal, nrTotDecolonized;
  /*
  ** As yet unused count variables.
  **
  ** int nrTotVegResRecover, nrTotSeedBankRecover, nrStepVegResilient,
  **     nrStepSeedBank, nrVegResilient, nrStepSeedBankRecover,
  **     nrStepVegResRecover, nrSeedBank, nrDecolonized;
  */

  /*
  ** Matrices:
  **   - currentState:   Values in [-32768;32767].
  **   - habSuitability: Values in [0;1000].
  **   - barriers:       Values in [0;255].
  **   - pixelAge:       Values in [0;255].
  **   - roadRiver:      Values in [0;255].
  **   - noDispersal:    Values in [0;255].
  */
  int **currentState, **habSuitability, **barriers, **pixelAge,
      **roadRiver, **noDispersal;

  /*
  ** Initialize the variables.
  */
  advOutput = false;
  currentState = NULL;
  habSuitability = NULL;
  barriers = NULL;
  pixelAge = NULL;
  roadRiver = NULL;
  noDispersal = NULL;
  seedProdProb = NULL;
  dispKernel = NULL;
  if (mcInit (*paramFile) == -1)
  {
    *nrFiles = -1;
    goto End_of_Routine;
  }
  startTime = time (NULL);
  /*
  ** We'll use 'srand' and 'rand' here, as 'srandom' and 'random' do not work
  ** on Windows :-(
  */
  srand (startTime);
  
  /*
  ** Load and prepare the data.
  **
  ** Species initial distribution.
  */
  currentState = (int **)malloc (nrRows * sizeof (int *));
  for (i = 0; i < nrRows; i++)
  {
    currentState[i] = (int *)malloc (nrCols * sizeof (int));
  }
  sprintf (fileName, "%s.asc", initDistrFile);
  if (mcReadMatrix (fileName, currentState) == -1)
  {
    *nrFiles = -1;
    goto End_of_Routine;
  }
  /*
  ** Habitat suitability matrix.
  */
  habSuitability = (int **)malloc (nrRows * sizeof (int *));
  for (i = 0; i < nrRows; i++)
  {
    habSuitability[i] = (int *)malloc (nrCols * sizeof (int));
  }
  /*
  ** Barrier options.
  */
  barriers = (int **)malloc (nrRows * sizeof (int *));
  for (i = 0; i < nrRows; i++)
  {
    barriers[i] = (int *)malloc (nrCols * sizeof (int));
    for (j = 0; j < nrCols; j++)
    {
      barriers[i][j] = 0;
    }
  }
  if (useBarrier)
  {
    sprintf (fileName, "%s.asc", barrierFile);
    if (mcReadMatrix (fileName, barriers) == -1)
    {
      *nrFiles = -1;
      goto End_of_Routine;
    }
    mcFilterByBarrier (currentState, barriers);
  }
  /* Time to reach dispersal maturity.
  **
  ** Before a pixel (= a population) can disperse, it needs to reach a certain
  ** user-defined age. The "age" is a matrix that keeps track of the age
  ** of all pixels:
  **   0 = pixel is not colonized.
  **   1 = pixel is colonized since 1 "dispersal event".
  **   2 = pixel is colonized since 2 "dispersal event".
  **   3 = etc...
  ** Fill the "age" to reflect the initial distribution of the species:
  ** where the species is present, pixels get a value of 'FullMaturity', where
  ** the species is absent, pixels get a value of 0.
  */
  pixelAge = (int **)malloc (nrRows * sizeof (int *));
  for (i = 0; i < nrRows; i++)
  {
    pixelAge[i] = (int *)malloc (nrCols * sizeof (int));
    for (j = 0; j < nrCols; j++)
    {
      if (currentState[i][j] == 1)
      {
    pixelAge[i][j] = fullMatAge;
      }
      else
      {
    pixelAge[i][j] = 0;
      }
    }
  }
  /*
  ** The "no dispersal" matrix.
  **
  ** This Matrix will keep track of the species distribution under the
  ** "no dispersal" scenario.
  */
  noDispersal = (int **)malloc (nrRows * sizeof (int *));
  for (i = 0; i < nrRows; i++)
  {
    noDispersal[i] = (int *)malloc (nrCols * sizeof (int));
    for (j = 0; j < nrCols; j++)
    {
      noDispersal[i][j] = currentState[i][j];
    }
  }
            
  /*
  ** Initialize counter variables.
  **
  ** Reset pixel counters to zero before we start the dispersal
  ** simulation.
  */
  nrInitial = 0;
  nrColonized = 0;
  nrAbsent = 0;
  nrNoDispersal = 0;
  nrUnivDispersal = 0;
  nrTotColonized = 0;
  nrTotDecolonized = 0;
  nrTotLDDSuccess = 0;
  tempResilience = true;
  nrStepColonized = 0;
  nrStepDecolonized = 0;
  nrStepLDDSuccess = 0;
  /*
  ** Currently unused.
  **
  ** nrVegResilient = 0;
  ** nrSeedBank = 0;
  ** nrTotVegResRecover = 0;
  ** nrTotSeedBankRecover = 0;
  */

  /*
  ** Count the number of initially colonized pixels
  ** (i.e. initial species distribution).
  */
  for (i = 0; i < nrRows; i++)
  {
    for (j = 0; j < nrCols; j++)
    {
      if (currentState[i][j] == 1)
      {
        nrInitial++;
      }
    }
  }
  nrColonized = nrInitial;
  nrNoDispersal = nrInitial;
  nrUnivDispersal = nrInitial;
  nrAbsent = (nrRows * nrCols) - nrInitial;

  /*
  ** Write the initial state to the data file.
  */
  sprintf (fileName, "%s/%s_stats.txt", simulName, simulName);
  if ((fp = fopen (fileName, "w")) != NULL)
  {
    fprintf (fp, "envChgStep\tdispStep\tloopID\tnrUnivDispersal\tnrNoDispersal\tnrColonized\tnrAbsent\tnrStepColonized\tnrStepDecolonized\tnrStepLDDSuccess\n");
    fprintf (fp, "0\t0\t0\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", nrUnivDispersal,
	     nrNoDispersal, nrColonized, nrAbsent, nrStepColonized,
	     nrStepDecolonized, nrStepLDDSuccess);
  }
  else
  {
    *nrFiles = -1;
    Rprintf ("Could not open statistics file for writing.\n");
    goto End_of_Routine;
  }
    
  /*
  ** Simulate plant dispersal and migration (the core of the method).
  **
  ** The environmental change loop (if simulation is run without change in
  ** environment this loop runs only once).
  */
  Rprintf ("Running MigClim.\n");
  for (envChgStep = 1; envChgStep <= nrEnvChgSteps; envChgStep++)
  {
    /*
    ** Print the current environmental change iteration.
    */
    Rprintf ("  %d...\n", envChgStep);

    /*
    ** Update and reclassify the habitat suitability matrix.
    ** Reclassify habSuitability according to the user-defined threshold.
    ** If the user selected to use "habitat invasibility" we do not reclassify
    ** values >= Threshold. If user selected not to use "habitat invasibility"
    ** then we rclassify values >= Threshold to 1000.
    */
    sprintf (fileName, "%s%d.asc", hsMapFile, envChgStep);
    if (mcReadMatrix (fileName, habSuitability) == -1)
    {
      *nrFiles = -1;
      goto End_of_Routine;
    }
    if (rcThreshold > 0)
    {
      for (i = 0; i < nrRows; i++)
      {
	for (j = 0; j < nrCols; j++)
	{
	  if (habSuitability[i][j] < rcThreshold)
	  {
	    habSuitability[i][j] = 0;
	  }
	  else
	  {
	    habSuitability[i][j] = 1000;
	  }
	}
      }
    }
    
    /*
    ** Set the values that will keep track of pixels colonized during the next
    ** climate change loop.
    ** "loopID" is the value that will be given to the pixel colonized
    ** during the current loop. The pixels being decolonized during the
    ** current loop will be given a value of "-loopID". The limits in
    ** terms of number of dispersal events and environmental change
    ** events is the following:
    **  - Without environmental change, the maximum number of dispStep
    **    is 29'500.
    **  - With environmental change, the maximum number of dispStep is
    **    98, the number of environmental change loops is limited to 250.
    ** The coding for the loopID is as follows: envChgStep * 100 + dispStep.
    */
    if (envChgStep == 0)
    {
      loopID = 1;
    }
    else
    {
      loopID = envChgStep * 100;
    }
      
    /*
    ** "Unlimited" and "no dispersal" scenario pixel count.
    ** Here we compute the number of pixels that would be colonized if dispersal
    ** was unlimited or null. This is simply the sum of all potentially suitable
    ** habitats (excluding barriers and filters).
    */
    nrUnivDispersal = mcUnivDispCnt (habSuitability, barriers);
    updateNoDispMat (habSuitability, noDispersal, &nrNoDispersal);

    /*
    ** Update for temporarily resilient pixels.
    */
    for (i = 0; i < nrRows; i++)
    {
      for (j = 0; j < nrCols; j++)
      {
	/*
	** Udate non-suitable pixels.
	** If a pixel turned unsuitable, we update its status to
	** "Temporarily Resilient".
	*/
	if ((habSuitability[i][j] == 0) && (currentState[i][j] > 0))
	{
	  nrStepDecolonized++;
	  /*
	  ** If the user selected TemporaryResilience==T, then the pixel is
	  ** set to "Temporary Resilient" status.
	  */
	  if (tempResilience == true)
	  {
	    currentState[i][j] = 29900;
	  }
	  else
	  {
	    /*
	    ** If not temporary resilience was specified, then the pixel is
	    ** set to "decolonized" status.
	    */
	    currentState[i][j] = -1 - loopID;
	    pixelAge[i][j] = 0;
	    /*
	    ** NOTE: Later we can add "Vegetative" and "SeedBank" resilience
	    ** options at this location.
	    */
	  }
	}
      }
    }

    /*
    **Dispersal event loop starts here.
    */
    for (dispStep = 1; dispStep <= nrDispSteps; dispStep++)
    {
      /*
      ** Set the value of "loopID" for the current iteration of the dispersal
      ** loop.
      */
      loopID++;
      
      /*
      ** Reset pixel counters that count pixels within the current loop.
      */
      nrStepColonized = 0;
      nrStepLDDSuccess = 0;
      if (dispStep > 1)
      {
          nrStepDecolonized = 0;
      }
      /*
      ** Currently unused variables.
      **
      ** nrStepVegResilient = 0;
      ** nrStepSeedBank = 0;
      ** nrStepVegResRecover = 0;
      ** nrStepSeedBankRecover = 0;
      */
      
      /*
      ** Source cell search: Can the sink pixel be colonized?
      ** There are 4 conditions to be met for a sink pixel to become colonized:
      **   1. Sink pixel is currently suitable and not already colonised.
      **   2. Sink pixel is within dispersal distance of an already colonised
      **      and mature pixel.
      **   3. Source pixel has reached dispersal maturity.
      **   4. There is no obstacle (barrier) between the pixel to be colonised
      **      (sink pixel) and the pixel that is already colonised (source
      **      pixel).
      **
      ** Loop through the cellular automaton.
      */
      for (i = 0; i < nrRows; i++)
      {
	for (j = 0; j < nrCols; j++)
	{
	  /*
	  ** Reset variables.
	  */
	  habIsSuitable = false;
	  cellInDispDist = false;
	  /*
	  ** 1. Test whether the pixel is a suitable sink (i.e., its habitat is
	  **    suitable, it's unoccupied and is not on a barrier or filter
	  **    pixel).
	  */
	  if ((habSuitability[i][j] > 0) && (currentState[i][j] <= 0) &&
	      (barriers[i][j] == 0))
	  {
	    habIsSuitable = true;
	  }
	  /*
	  ** 2. Test whether there is a source cell within the dispersal
	  **    distance. To be more time efficient, this code runs only if the
	  **    answer to the first question is positive. Additionally, if there
	  **    is a source cell within dispersion distance and if a barrier was
	  **    asked for, then we also check that there is no barrier between
	  **    this source cell and the sink cell (this verification is carried
	  **    out in the "SearchSourceCell" function).
	  */
	  if (habIsSuitable)
	  {
	    /*
	    ** Now we search if there is a suitable source cell to colonize the
	    ** sink cell.
	    */
	    if (mcSrcCell (i, j, currentState, pixelAge, loopID,
	        habSuitability[i][j], barriers))
	    {
	        cellInDispDist = true;
	    }
	  }
	  /*
	  ** Update pixel status.
	  */
	  if (habIsSuitable && cellInDispDist)
	  {
	    /*
	    ** Only if the 2 conditions are fullfilled the cell's status is set
	    ** to colonised.
	    */
	    currentState[i][j] = loopID;
	    nrStepColonized++;
	    /* 
	    ** If the pixel was in seed bank resilience state, then we
	    ** update the corresponding counter. Currently not used.
	    **
	    ** if (pixelAge[i][j] == 255)
	    ** {
	    **   nrStepSeedBank--;
	    ** }
	    */
	      
	    /*
	    ** Update "age" value. We do this only now because we needed the
	    ** old "age" value just before to determine whether a pixel was in
	    ** "Decolonized" or "SeedBank resilience" status.
	    */
	    pixelAge[i][j] = 0;
	  }
	}
      }

      /*
      ** If the LDD frequence is larger than zero, perform it.
      */
      if (lddFreq > 0.0)
      {
	/*
	** Loop through the entire cellular automaton.
	*/
	for (i = 0; i < nrRows; i++)
	{
	  for (j = 0; j < nrCols; j++)
	  {
	    /*
	    ** Check if the pixel is a source cell (i.e. is it colonised since
	    ** at least 1 dispersal Loop) and check if the pixel has reached
	    ** dispersal maturity.
	    */
	    if ((currentState[i][j]) > 0 && (currentState[i][j] != loopID))
	    {
	      if (pixelAge[i][j] >= initMatAge)
	      {
		/*
		** Set the probability of generating an LDD event. This
		** probability is weighted by the age of the cell.
		*/
		if (pixelAge[i][j] >= fullMatAge)
		{
		  lddSeedProb = lddFreq;
		}
		else
		{
		  lddSeedProb = lddFreq * seedProdProb[pixelAge[i][j] -
		                initMatAge];
		}
		/*
		** Now we can try to generate a LDD event with the calculated
		** probability.
		*/
		if (UNIF01 < lddSeedProb)
		{
		  /*
		  ** Randomly select a pixel within the distance
		  ** "minDist - maxDist".
		  */
		  mcRandomPixel (&rndPixel);
		  rndPixel.row = rndPixel.row + i;
		  rndPixel.col = rndPixel.col + j;
		  /*
		  ** Now we check if this random cell is a suitable sink cell.
		  */
		  if (mcSinkCellCheck (rndPixel, currentState, habSuitability,
		                      barriers))
		  {
		    /*
		    ** The pixel gets colonized.
		    */
		    currentState[rndPixel.row][rndPixel.col] = loopID;
		    nrStepColonized++;
		    nrStepLDDSuccess++;
		    /* 
		    ** If the pixel was in seed bank resilience state, then we
		    ** update the corresponding counter. Currently not used.
		    **
		    ** if (pixelAge[rndPixel.row][rndPixel.col] == 255)
		    ** {
		    **   nrStepSeedBank--;
		    ** }
		    */
		    
		    /*
		    ** Reset pixel age.
		    */
		    pixelAge[rndPixel.row][rndPixel.col] = 0;                  
		  }
		}
	      }
	    }
	  }
	}
      }
            
      /*
      ** Update pixel age:
      ** At the end of a dispersal loop we want to
      ** increase the "age" of each colonized pixel.
      **
      ** Reminder: pixel "age" structure is as follows:
      **   0 = Pixel is either "Absent", "Decolonized" or has just been
      **       "Colonized" during this dispersal step.
      **   1 to 250 = Pixel is in "Colonized" or "Temporarily Resilient" status.
      **       The value indicates the number of "dispersal events (usually
      **       years) since when the pixel was colonized.
      **   255 = Pixel is in "SeedBank Resilience" state.
      */
      for (i = 0; i < nrRows; i++)
      {
	for (j = 0; j < nrCols; j++)
	{
	  /*
	  ** If the pixel is in "Colonized" or "Temporarily Resilient" state,
	  ** update it's age value.
	  */
	  if (currentState[i][j] > 0)
	  {
	    pixelAge[i][j] += 1;
	  }
	  /*
	  ** If a pixel is in "Temporarily Resilient" state, we also increase
	  ** its "currentState" value by 1, so that the pixels gains 1
	  ** year of "Temporarily Resilience" age.
	  */
	  if (currentState[i][j] >= 29900)
	  {
	    currentState[i][j] += 1;
	  }
	}
      }

      /*
      ** Update pixel counters.
      */
      nrColonized = nrColonized + nrStepColonized - nrStepDecolonized;
      nrAbsent = (nrRows * nrCols) - nrColonized;
      nrTotColonized += nrStepColonized;
      nrTotDecolonized += nrStepDecolonized;
      nrTotLDDSuccess += nrStepLDDSuccess;
      /*
      ** Currently unused variables.
      **
      ** nrTotVegResRecover += nrStepVegResRecover;
      ** nrTotSeedBankRecover += nrStepSeedBankRecover;
      */
      
      /*
      ** Write current iteration data to the statistics file.
      */
      fprintf (fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", envChgStep,
           dispStep, loopID, nrUnivDispersal, nrNoDispersal, nrColonized,
           nrAbsent, nrStepColonized, nrStepDecolonized, nrStepLDDSuccess);

      /*
      ** If the user has requested full output, also write the current state
      ** matrix to file.
      */
      if (fullOutput)
      {
	sprintf (fileName, "%s/%s_step_%d.asc", simulName, simulName, loopID);
	if (mcWriteMatrix (fileName, currentState) == -1)
	{
	  *nrFiles = -1;
	  goto End_of_Routine;
	}
      }
    } /* END OF: dispStep */
    
    /*
    ** Update temporarily resilient pixels.
    ** Temporarily resilient pixels can be distinguished by:
    **   - CurrentState_Matrix = 29'900 to 29'999. Increases by 1 at each
    **     year. Age_Matrix has a positive value.
    */
    for (i = 0; i < nrRows; i++)
    {
      for (j = 0; j < nrCols; j++)
      {
	if (currentState[i][j] >= 29900)
	{
	  currentState[i][j] = nrDispSteps - loopID - 1;
	  pixelAge[i][j] = 0;
	}
      }
    }
  } /* END OF: envChgStep */
  Rprintf ("done.\n");
  
  /*
  ** Update currentState matrix for pixels that are suitable but
  ** could not be colonized due to dispersal limitations.
  ** These pixels are assigned a value of 30'000
  */
  for (i = 0; i < nrRows; i++)
  {
    for (j = 0; j < nrCols; j++)
    {
      if ((habSuitability[i][j] > 0) && (currentState[i][j] <= 0) &&
	  (barriers[i][j] == 0))
      {
	currentState[i][j] = 30000;
      }
    }
  }
  
  /*
  ** Write the final state matrix to file.
  */
  sprintf (fileName, "%s/%s_raster.asc", simulName, simulName);
  if (mcWriteMatrix (fileName, currentState) == -1)
  {
    *nrFiles = -1;
    goto End_of_Routine;
  }
  
  /*
  ** Write summary output to file.
  */
  simulTime = time (NULL) - startTime;
  sprintf (fileName, "%s/%s_summary.txt", simulName, simulName);
  if ((fp2 = fopen (fileName, "w")) != NULL)
  {
    fprintf (fp2, "SimulationName\tIniCount\tNoDispersalCount\tUnlimitedDispersalCount\tDispersalCount\tAbsentCount\tnrTotColonized\tnrTotDecolonized\tnrTotLDDSuccess\tRunTime\n");
    fprintf (fp2, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", simulName,
	     nrInitial, nrNoDispersal, nrUnivDispersal, nrColonized, nrAbsent,
	     nrTotColonized, nrTotDecolonized, nrTotLDDSuccess, simulTime);
    fclose (fp2);
  }
  else
  {
    *nrFiles = -1;
    Rprintf ("Could not write summary output to file.\n");
    goto End_of_Routine;
  }
  
  /*
  ** Set the number of output files created.
  */
  *nrFiles = nrEnvChgSteps;

 End_of_Routine:
  /*
  ** Close the data file.
  */
  if (fp != NULL)
  {
    fclose (fp);
  }
  /*
  ** Free the allocated memory.
  */
  if (currentState != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (currentState[i]);
    }
    free (currentState);
  }
  if (habSuitability != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (habSuitability[i]);
    }
    free (habSuitability);
  }
  if (barriers != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (barriers[i]);
    }
    free (barriers);
  }
  if (pixelAge != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (pixelAge[i]);
    }
    free (pixelAge);
  }
  if (roadRiver != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (roadRiver[i]);
    }
    free (roadRiver);
  }
  if (noDispersal != NULL)
  {
    for (i = 0; i < nrRows; i++)
    {
      free (noDispersal[i]);
    }
    free (noDispersal);
  }
  if (dispKernel != NULL)
  {
    free (dispKernel);
  }
  if (seedProdProb != NULL)
  {
    free (seedProdProb);
  }
  /*
  ** Check the final result.
  */
  if (*nrFiles == -1)
  {
    Rprintf ("MigClim simulation aborted...\n");
  }
}    


/*
** mcRandomPixel: Select a random pixel from a central point (0;0) and within a
**                radius of at least minDist and at most maxDist.
**
** Parameters:
**   - pix:     A pointer to the pixel data structure to put the row and col in.
*/

void mcRandomPixel (pixel *pix)
{
  double rndDist, rndAngle;
  
  /*
  ** Select a random distance between minDist and maxDist, and a random angle
  ** between 0 and 2*pi.
  */
  rndDist = (UNIF01 * (maxDist - minDist)) + minDist;
  rndAngle = UNIF01 * 6.283185;
    
  /*
  ** Convert distance and angle into pixel row and column values.
  */
  pix->row = (int)(rndDist * cos (rndAngle));
  pix->col = (int)(rndDist * sin (rndAngle));
}


/*
** mcSinkCellCheck: Perform a basic check to see whether a given pixel fulfills
**                  the conditions to be a "sink" pixel (i.e. a pixel to be
**                  potentially colonized).
**
** The conditions that are checked are the following:
**   1° Pixel must be within the limits of the cellular automaton.
**   2° Pixel must be empty (i.e. non-occupied).
**   3° Pixel must contain suitable habitat (i.e. it must be potentially
**      suitable).
**   4° Pixel does not belong to a barrier or filter.
**
** Parameters:
**   - pix: The pixel to consider.
**   - curState: A pointer to the current state matrix.
**   - habSuit:  A pointer to the habitat suitability matrix.
**   - bars:     A pointer to the barriers matrix.
**
** Returns:
**   If the pixel is suitable: true.
**   Otherwise:                false.
*/

bool mcSinkCellCheck (pixel pix, int **curState, int **habSuit, int **bars)
{
  bool suitable;

  suitable = true;

  /*
  ** 1° Verify the pixel is within the limits of the cellular automaton.
  */
  if ((pix.row < 0) || (pix.row >= nrRows) || (pix.col < 0) ||
      (pix.col >= nrCols))
  {
    suitable = false;
  }
  /*
  ** 2° Verify the pixel is empty.
  */
  else if (curState[pix.row][pix.col] > 0)
  {
    suitable = false;
  }
  /*
  ** 3° Verify the pixel contains suitable habitat for the species.
  */
  else if ((UNIF01 * 1000) > habSuit[pix.row][pix.col])
  {
    suitable = false;
  }
  /*
  ** 4° Verify the pixel is not a "barrier" pixel.
  */
  else if (bars[pix.row][pix.col] > 0)
  {
    suitable = false;
  }
    
  /*
  ** Return the result.
  */
  return (suitable);
}


/*
** EoF: migrate.c
*/
