/*
** file_io.c: Functions for file input/output.
**
** Wim Hordijk   Last modified: 05 October 2011
*/

#include "migclim.h"


/*
** mcInit: Initialize the MigClim model by reading the parameter values from
**         file.
**
** Parameters:
**   - paramFile: The name of the file from which to read the parameter values.
**
** Returns:
**   - If everything went fine:  0.
**   - Otherwise:               -1.
*/

int mcInit (char *paramFile)
{
  int    i, age, status, lineNr;
  char   line[1024], param[64];
  float  p;
  FILE  *fp;

  status = 0;
  fp = NULL;

  /*
  ** Set default parameter values.
  */
  nrRows = 0;
  nrCols = 0;
  strcpy (initDistrFile, "");
  strcpy (hsMapFile, "");
  strcpy (barrierFile, "");
  useBarrier = false;
  barrierType = WEAK_BARRIER;
  nrEnvChgSteps = 0;
  nrDispSteps = 0;
  dispDist = 0;
  initMatAge = 0;
  fullMatAge = 0;
  rcThreshold = 0;
  minDist = 10;
  maxDist = 15;
  lddFreq = 0.0;
  fullOutput = false;
  strcpy (simulName, "mySimul");
  
  /*
  ** Open the file for reading.
  */
  if ((fp = fopen(paramFile, "r")) == NULL)
  {
    status = -1;
    Rprintf ("Can't open parameter file %s\n", paramFile);
    goto End_of_Routine;
  }

  /*
  ** While there are lines left, read and parse them.
  */
  lineNr = 0;
  param[0] = '\0';
  while (fgets (line, 1024, fp) != NULL)
  {
    lineNr++;
    sscanf (line, "%s", param);
    /* nrRows */
    if (strcmp (param, "nrRows") == 0)
    {
      if ((sscanf (line, "nrRows %d", &nrRows) != 1) || (nrRows < 1))
      {
	status = -1;
	Rprintf ("Invalid number of rows on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* nrCols */
    else if (strcmp (param, "nrCols") == 0)
    {
      if ((sscanf (line, "nrCols %d", &nrCols) != 1) || (nrCols < 1))
      {
	status = -1;
	Rprintf ("Invalid number of columns on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* initDistrFile */
    else if (strcmp (param, "initDistrFile") == 0)
    {
      if (sscanf (line, "initDistrFile %s", initDistrFile) != 1)
      {
	status = -1;
	Rprintf ("Invalid initial distribution file name on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* hsMapFile */
    else if (strcmp (param, "hsMapFile") == 0)
    {
      if (sscanf (line, "hsMapFile %s", hsMapFile) != 1)
      {
	status = -1;
	Rprintf ("Invalid habitat suitability map file name on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* barrierFile */
    else if (strcmp (param, "barrierFile") == 0)
    {
      if (sscanf (line, "barrierFile %s", barrierFile) != 1)
      {
	status = -1;
	Rprintf ("Invalid barrier file name on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
      useBarrier = true;
    }
    /* barrierType */
    else if (strcmp (param, "barrierType") == 0)
    {
      if (sscanf (line, "barrierType %s", param) != 1)
      {
	status = -1;
	Rprintf ("Invalid barrier type on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
      if (strcmp (param, "weak") == 0)
      {
	barrierType = WEAK_BARRIER;
      }
      else if (strcmp (param, "strong") == 0)
      {
	barrierType = STRONG_BARRIER;
      }
      else
      {
	status = -1;
	Rprintf ("Invalid barrier type on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* nrEnvChgSteps */
    else if (strcmp (param, "nrEnvChgSteps") == 0)
    {
      if ((sscanf (line, "nrEnvChgSteps %d", &nrEnvChgSteps) != 1) ||
	  (nrEnvChgSteps < 1))
      {
	status = -1;
	Rprintf ("Invalid number of environmental change steps on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* nrDispSteps */
    else if (strcmp (param, "nrDispSteps") == 0)
    {
      if ((sscanf (line, "nrDispSteps %d", &nrDispSteps) != 1) ||
	  (nrDispSteps < 1))
      {
	status = -1;
	Rprintf ("Invalid number of dispersal steps on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* dispDist & dispKernel */
    else if (strcmp (param, "dispDist") == 0)
    {
      if ((sscanf (line, "dispDist %d", &dispDist) != 1) ||
	  (dispDist < 1))
      {
	status = -1;
	Rprintf ("Invalid dispersal distance on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
      lineNr++;
      dispKernel = (double *)malloc (dispDist * sizeof (double));
      if (fscanf (fp, "dispKernel %f", &p) != 1)
      {
	status = -1;
	Rprintf ("Dispersal kernel expected on line %d in parameter file %s\n",
		 line, paramFile);
	goto End_of_Routine;
      }
      dispKernel[0] = p;
      for (i = 1; i < dispDist; i++)
      {
	if (fscanf (fp, "%f", &p) != 1)
	{
	  status = -1;
	  Rprintf ("Invalid dispersal kernel values on line %d in parameter file %s.\n",
		   line, paramFile);
	  goto End_of_Routine;
	}
	dispKernel[i] = p;
      }
      p = fscanf (fp, "\n");
    }
    /* initMatAge */
    else if (strcmp (param, "initMatAge") == 0)
    {
      if ((sscanf (line, "initMatAge %d", &initMatAge) != 1) ||
	  (initMatAge < 1))
      {
	status = -1;
	Rprintf ("Invalid initial maturity age on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* fullMatAge */
    else if (strcmp (param, "fullMatAge") == 0)
    {
      if ((sscanf (line, "fullMatAge %d", &fullMatAge) != 1) ||
	  (fullMatAge < 1))
      {
	status = -1;
	Rprintf ("Invalid full maturity age on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
      lineNr++;
      age = fullMatAge - initMatAge;
      if (age == 0)
      {
	age = 1;
      }
      seedProdProb = (double *)malloc (age * sizeof (double));
      if (fscanf (fp, "seedProdProb %f", &p) != 1)
      {
	status = -1;
	Rprintf ("Seed production probabilities expected on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
      seedProdProb[0] = p;
      for (i = 1; i < age; i++)
      {
	if (fscanf (fp, "%f", &p) != 1)
	{
	  status = -1;
	  Rprintf ("Invalid seed production probability on line %d in parameter file %s\n",
		   lineNr, paramFile);
	  goto End_of_Routine;
	}
	seedProdProb[i] = p;
      }
      p = fscanf (fp, "\n");
    }
    /* rcThreshold */
    else if (strcmp (param, "rcThreshold") == 0)
    {
      if ((sscanf (line, "rcThreshold %d", &rcThreshold) != 1) ||
	  (rcThreshold < 0) || (rcThreshold > 1000))
      {
	status = -1;
	Rprintf ("Invalid reclassification threshold on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* lddFreq */
    else if (strcmp (param, "lddFreq") == 0)
    {
      if ((sscanf (line, "lddFreq %f", &p) != 1) ||
	  (p < 0.0) || (p > 1.0))
      {
	status = -1;
	Rprintf ("Invalid long-distance dispersal frequency on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
      lddFreq = p;
    }
    /* minDist */
    else if (strcmp (param, "minDist") == 0)
    {
      if ((sscanf (line, "minDist %d", &minDist) != 1) || (minDist < 0))
      {
	status = -1;
	Rprintf ("Invalid minimum long-distance dispersal value on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* maxDist */
    else if (strcmp (param, "maxDist") == 0)
    {
      if ((sscanf (line, "maxDist %d", &maxDist) != 1) || (maxDist <= minDist))
      {
	status = -1;
	Rprintf ("Invalid maximum long-distance dispersal value on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* fullOutput */
    else if (strcmp (param, "fullOutput") == 0)
    {
      if (sscanf (line, "fullOutput %s", param) != 1)
      {
	status = -1;
	Rprintf ("Incomplete 'fullOutput' argument on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
      if (strcmp (param, "true") == 0)
      {
	fullOutput = true;
      }
      else if (strcmp (param, "false") == 0)
      {
	fullOutput = false;
      }
      else
      {
	status = -1;
	Rprintf ("Invalid value for argument 'fullOutput' on line %d in parameter file %s\n", lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* simulName */
    else if (strcmp (param, "simulName") == 0)
    {
      if (sscanf (line, "simulName %s", simulName) != 1)
      {
	status = -1;
	Rprintf ("Invalid output file name on line %d in parameter file %s\n",
		 lineNr, paramFile);
	goto End_of_Routine;
      }
    }
    /* Unknown parameter */
    else
    {
      status = -1;
      Rprintf ("Unknown parameter on line %d in parameter file %s\n", lineNr,
	      paramFile);
      goto End_of_Routine;
    }
  }

  /*
  ** Check the parameter values for validity.
  */
  if (nrRows == 0)
  {
    status = -1;
    Rprintf ("No number of rows specified in parameter file %s\n", paramFile);
    goto End_of_Routine;
  }
  if (nrCols == 0)
  {
    status = -1;
    Rprintf ("No number of columns specified in parameter file %s\n",
	     paramFile);
    goto End_of_Routine;
  }
  if (strlen (initDistrFile) == 0)
  {
    status = -1;
    Rprintf ("No initial distribution file name specified in parameter file %s\n",
	     paramFile);
    goto End_of_Routine;
  }
  if (strlen (hsMapFile) == 0)
  {
    status = -1;
    Rprintf ("No habitat suitability map file name specified in parameter file %s\n",
	     paramFile);
    goto End_of_Routine;
  }
  if (useBarrier && (strlen (barrierFile) == 0))
  {
    status = -1;
    Rprintf ("No barrier file name specified in parameter file %s\n",
	     paramFile);
    goto End_of_Routine;
  }
  if (nrEnvChgSteps == 0)
  {
    status = -1;
    Rprintf ("No number of environmental change steps specified in parameter file %s\n",
	     paramFile);
    goto End_of_Routine;
  }
  if (nrDispSteps == 0)
  {
    status = -1;
    Rprintf ("No number of dispersal steps specified in parameter file %s\n",
	     paramFile);
    goto End_of_Routine;
  }
  if (dispDist == 0)
  {
    status = -1;
    Rprintf ("No dispersal distance specified in parameter file %s\n",
	     paramFile);
    goto End_of_Routine;
  }
  if (initMatAge == 0)
  {
    status = -1;
    Rprintf ("No initial maturity age specified in parameter file %s\n",
	     paramFile);
    goto End_of_Routine;
  }
  if (fullMatAge == 0)
  {
    status = -1;
    Rprintf ("No full maturity age specified in parameter file %s\n",
	     paramFile);
    goto End_of_Routine;
  }
  
 End_of_Routine:
  /*
  ** Close the file and return the status.
  */
  if (fp != NULL)
  {
    fclose (fp);
  }
  return (status);
}


/*
** mcReadMatrix: Read a data matrix from file.
**
** Parameters:
**   - fname:  The name of the file to read from.
**   - mat:    The matrix to put the data in (assumed to be large enough).
**
** Returns:
**   - If everything went fine:  0.
**   - Otherwise:               -1.
*/

int mcReadMatrix (char *fname, int **mat)
{
  int   i, j, val, status, noData;
  char  line[1024], param[128];
  FILE *fp;

  status = 0;
  fp = NULL;
  
  /*
  ** Open the file for reading.
  */
  if ((fp = fopen(fname, "r")) == NULL)
  {
    status = -1;
    Rprintf ("Can't open data file %s\n", fname);
    goto End_of_Routine;
  }

  /*
  ** Get the 'meta data'.
  */
  if ((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %d\n", param, &val) != 2) ||
      (strcasecmp (param, "ncols") != 0))
  {
    status = -1;
    Rprintf ("'ncols' expected in data file %s.\n", fname);
    goto End_of_Routine;
  }
  if (val != nrCols)
  {
    status = -1;
    Rprintf ("Invalid number of columns in data file %s\n", fname);
    goto End_of_Routine;
  }
  if ((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %d\n", param, &val) != 2) ||
      (strcasecmp (param, "nrows") != 0))
  {
    status = -1;
    Rprintf ("'nrows' expected in data file %s\n", fname);
    goto End_of_Routine;
  }
  if (val != nrRows)
  {
    status = -1;
    Rprintf ("Invalid number of rows in data file %s.\n", fname);
    goto End_of_Routine;
  }
  if ((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %d\n", param, &val) != 2) ||
      (strcasecmp (param, "xllcorner") != 0))
  {
    status = -1;
    Rprintf ("'xllcorner' expected in data file %s\n", fname);
    goto End_of_Routine;
  }
  if ((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %d\n", param, &val) != 2) ||
      (strcasecmp (param, "yllcorner") != 0))
  {
    status = -1;
    Rprintf ("'yllcorner' expected in data file %s\n", fname);
    goto End_of_Routine;
  }
  if ((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %d\n", param, &val) != 2) ||
      (strcasecmp (param, "cellsize") != 0))
  {
    status = -1;
    Rprintf ("'cellsize' expected in data file %s\n", fname);
    goto End_of_Routine;
  }
  if ((fgets (line, 1024, fp) == NULL) ||
      (sscanf (line, "%s %d\n", param, &noData) != 2) ||
      (strcasecmp (param, "nodata_value") != 0))
  {
    status = -1;
    Rprintf ("'nodata_value' expected in data file %s\n", fname);
    goto End_of_Routine;
  }
  
  /*
  ** Read the values into the matrix.
  */
  for (i = 0; i < nrRows; i++)
  {
    for (j = 0; j < nrCols; j++)
    {
      if (fscanf(fp, "%d", &val) != 1)
      {
	status = -1;
	Rprintf ("Invalid value in data file %s\n", fname);
	goto End_of_Routine;
      }
      if (val != noData)
      {
	mat[i][j] = val;
      }
      else
      {
	mat[i][j] = 0;
      }
    }
    val = fscanf (fp, "\n");
  }

  /*
  ** Close the file and return the status.
  */
 End_of_Routine:
  if (fp != NULL)
  {
    fclose (fp);
  }
  return (status);
}


/*
** mcWriteMatrix: Write a data matrix to file.
**
** Parameters:
**   - fname:  The name of the file to write to.
**   - mat:    The data matrix to write.
**
** Returns:
**   - If everything went fine:  0.
**   - Otherwise:               -1.
*/

int mcWriteMatrix (char *fname, int **mat)
{
  int   i, j, status;
  FILE *fp;

  status = 0;
  fp = NULL;
  
  /*
  ** Open the file for writing.
  */
  if ((fp = fopen(fname, "w")) == NULL)
  {
    status = -1;
    Rprintf ("Can't open data file %s for writing.\n", fname);
    goto End_of_Routine;
  }

  /*
  ** Write the 'meta data'.
  */
  fprintf (fp, "ncols %d\n", nrCols);
  fprintf (fp, "nrows %d\n", nrRows);
  fprintf (fp, "xllcorner 0\n");
  fprintf (fp, "yllcorner 0\n");
  fprintf (fp, "cellsize 10\n");
  fprintf (fp, "nodata_value -9999\n");
  
  /*
  ** Write the data to file.
  */
  for (i = 0; i < nrRows; i++)
  {
    for (j = 0; j < nrCols; j++)
    {
      fprintf(fp, "%d ", mat[i][j]);
    }
    fprintf (fp, "\n");
  }

  /*
  ** Close the file and return the status.
  */
 End_of_Routine:
  if (fp != NULL)
  {
    fclose (fp);
  }
  return (status);
}


/*
** EoF: file_io.c
*/
