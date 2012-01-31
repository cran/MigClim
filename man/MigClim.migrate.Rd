\name{MigClim.migrate}
\alias{MigClim.migrate}
\title{Initialize and run the MigClim simulation}
\description{Initialize the main MigClim simulation by passing the relevant
  parameter values, which will be written to a parameter file. The
  simulation is then run with these parameter values.}
\usage{MigClim.migrate (initDistrFile="InitialDist", hsMapFile="HSmap",
  barrierFile="", barrierType="weak", nrEnvChgSteps=1, nrDispSteps=1,
  dispKernel=c(1.0,1.0), initMatAge=1, fullMatAge=1, seedProdProb=c(1.0),
  rcThreshold=0, lddFreq=0.0, minDist=NULL, maxDist=NULL, fullOutput=FALSE,
  simulName="mySimul", overWrite=FALSE, testMode=FALSE, keepTempFiles=FALSE)}
\arguments{
  \item{initDistrFile}{The initial distribution of the species. This can be given either a string indicating the name of a raster object (see 'Details' for supported formats) or as a data frame object (see 'Details' for how to structure your data frame). Please note that the inputs for 'initDistrFile', 'hsMapFile' and 'barrierFile' (optional) must always be given in the same format.}
  \item{hsMapFile}{The habitat suitability values. This can be given as a string indicating the 'base' name of the raster files that contain the habitat suitability maps. Iteration numbers (1,2,3,...) are automatically added to this 'base' name to get the file name for the habitat suitability map for each next environmental change iteration (see the 'Details' section for supported formats). Alternatively, the habitat suitability information can also be given as a data frame object, where each column indicates a successive habitat suitability map (see the 'Details' section for further information on how this data frame must be structured).}
  \item{barrierFile}{The name of the raster file that contains barrier information or a single column data frame (or vector) containing this information. If an empty string is given (default value), no barrier information is used.}
  \item{barrierType}{The barrier type to use. Values can be either 'weak' (default value) or 'strong'. Not relevant if barrier information is not used. 'weak' barriers will allow dispersal to proceed through two diagonally adjacent barrier pixels, 'strong' barriers won't.}
  \item{nrEnvChgSteps}{The number of environmental change steps to perform. At each environmental change step the habitat suitability values are updated with the values of the corresponding habitat suitability map (and therefore the number of environmental change steps must match the number of habitat suitability maps available).}
  \item{nrDispSteps}{The number of dispersal steps to perform within each environmental change step. For instance, if one wants to simulate dispersal to occur once a year, and the habitat suitability maps represent 5 years intervals, then 'nrDispSteps' should be set to 5.}
  \item{dispKernel}{The dispersal kernel. A vector of dispersal probabilities (values between 0.0 and 1.0) giving the probability for a source cell to colonize an empty cell. The distance unit is the 'pixel', with the first value in the vector representing the probability for source cell to colonize a directly adjacent cell.}
  \item{initMatAge}{The initial maturity age. Newly colonized cells younger than this age cannot colonize other cells. The time unit is a dispersal step (which usually should be equal to 1 year).}
  \item{fullMatAge}{The full maturity age. Cells that have reached their full maturity can produce propagules with a probability of 1.}
  \item{seedProdProb}{The seed production probability. A vector of length 'fullMatAge'-'initMatAge' (or length 1 if both ages are equal) with seed production probabilities (values between 0.0 and 1.0).}
  \item{rcThreshold}{The reclassification threshold (a value between 0 and 1000; default=0). If 'rcThreshold > 0', then continuous values of the habitat suitability maps will be reclassified according to 'rcThreshold'. Values of habitat suitability < 'rcThreshold' are reclassifier to '0' (unsuitable habitat) and values >= 'rcThreshold' are reclassified to '1000' (fully suitable habitat).  In the case where 'rcThreshold=0', the habitat suitability values are not reclassified, and are instead considered as representing habitat 'invasibility', modulating the probability of an unoccupied cell to become colonized (suitability values are considered as probabilities * 1000).}
  \item{lddFreq}{The long-distance dispersal frequency, i.e., the probability of a long-distance dispersal event per (suitable) pixel. If set to 0.0 (default), no long-distance dispersal is performed. Value should be between 0.0 and 1.0.}
  \item{minDist}{The minimum distance for long-distance dispersal (not used if 'lddFreq' is 0.0). This value must be larger than the regular dispersal ditance(given by the length of 'dispKernel').}
  \item{maxDist}{The maximum distance for long-distance dispersal (not used if 'lddFreq' is 0.0).}
  \item{fullOutput}{If 'TRUE', the current state of the simulation is written to am ACII raster file after each dispersal step (allowing to reconstruct the dispersal process at each step). If 'FALSE' (default), only the final state of the simulation is written to an ASCII grid file.}
  \item{simulName}{The 'base' name of for the different outputs produced by the MigClim simulation. Three different outputs are produced by the MigClim simulation (i.e, after running the 'MigClim.migrate' function): an ESRI ASCII Raster file named 'simulName'+'_raster.asc' that contains the final state of the simulation, an 'simulName'+'_stats.txt' file that contains the simulation's outputs after each dispersal event, and an 'simulName'+'_summary.txt' file that contains a single-line summary of the entire simulation.}
  \item{overWrite}{If 'TRUE' than any existing file with the same name as the ouputs of MigClim.migrate will be overwritten. If 'FALSE' then the function will stop if any ouput file does already exist.}
  \item{testMode}{If TestMode is set to 'TRUE' then the MigClim.migrate function will check all the provided input data but will not run the actual simulation. Use this for testing your data before running simulations that might take a long time or that are part of a loop.}
  \item{keepTempFiles}{If 'FALSE' (default), then any ".asc" file created from a conversion process in the function will be deleted when simulation is completed. If you wish to keep these files then set the value of this parameter to 'TRUE'.}
}
\details{The input data for initial distribution ('initDistrFile'), habitat suitability ('hsMapFile'), and (optionally) barriers ('barrierFile') can be provided as either a string giving the name of a raster file (located in the working directory) or as a data frame object. For a given simulation, all these inputs must be given in the same format.
Option 1: Input data provided as raster file. In this case, the input must be a string that contains the name of the raster file that is located in the working directory. The following raster formats are supported: (i) ascii GRID (files must have a '.asc' extension), (ii) R rasterLayer (see 'raster' package), (iii) ESRI GRID, (iv) GeoTIFF (files must have a '.tif' extension). Note that all input grids need to have exactly the same size (i.e., the same number of rows and the same number of columns).

The standard ESRI ASCII Raster format looks as follows (actual values depend on file content):
\preformatted{ncols         100
nrows         100
xllcorner     0
yllcorner     0
cellsize      10
nodata_value -9999
0 0 1 0 0 1 1 1 0 ...
1 0 1 0 1 0 1 0 1 ...
:
:
}

Option 2: Input data provided as data frames. It is important to note that, unlike when inputing data in a raster format, inputing it as data frames requires to provide the actual data frame object to the function, not a string indicating its name. 'initDistrFile' shoud have 3 colums (in this order): X coordinate, Y coordinate and initial distribution of the species (only values of 0 or 1 are accpted. 1 = species present, 0 = species absent.
\preformatted{Xcoordinate  Ycoordinate IniDist
520000       200100      1
520050       200100      1
520100       200100      0
...
}
          'hsMapFile' shoud have a number of colums equal to 'nrEnvChgSteps'. Each column giving the suitability of the habitat for a given environmental change step (in chronological order). Values of habitat suitability must be in the range [0:1000]
\preformatted{
HSmap1  HSmap2  HSmap3  ...
789     800     230     ...
560     720     150     ...
120     540      86     ...
...
}

'barrierFile' is optional and must have only one column (it can also be a vector). It must contain only values of either 1 (pixel is a barrier feature) or 0 (pixel is not a barrier feature).

The function output(s) will be written in ascii GRID format.
}
\value{The number of environmental change steps performed. The function also writes the following outputs into the current working directory: an ESRI ASCII Raster file named 'simulName'+'_raster.asc' that contains the final state of the simulation, an 'simulName'+'_stats.txt' file that contains the simulation's outputs after each dispersal event, and an 'simulName'+'_summary.txt' file that contains a single-line summary of the entire simulation.  state matrix is written as an ESRI ASCII Raster file with the file name given as the value of the 'outFile' argument in the MigClim.init function, with the extension '_raster.asc' added. The output file 'simulName' + '_stats.txt' contains summary statistics for each individual dispersal step. The output file 'simulName' + '_summary.txt' contains summary statistics over the entire simulation (including the running time in seconds).}
\references{Engler R., Guisan A., 2009. MIGCLIM: Predicting plant distribution and dispersal in a changing climate. Diversity and Distributions 15(4), 590-601.}
\seealso{MigClim.plot ()}
\examples{
  data(MigClim.TestData)

  ### Run MigClim with a data frame type input.
  n<-MigClim.migrate (initDistrFile=MigClim.TestData[,1:3],
     hsMapFile=MigClim.TestData[,4:8], barrierFile=MigClim.TestData[,9],
     barrierType="weak", nrEnvChgSteps=5, nrDispSteps=5,
     dispKernel=c(1.0,0.4,0.16,0.06,0.03), initMatAge=1, fullMatAge=5,
     seedProdProb=c(0.01,0.08,0.5,0.92), lddFreq=0.1, minDist=6,
     maxDist=15, fullOutput=FALSE, simulName="mcTest", 
     overWrite=TRUE, testMode=FALSE, keepTempFiles=FALSE)
  
\dontrun{  ### Convert data from a dataframe to ASCII grid format
  library(SDMTools)
  dataframe2asc(MigClim.TestData[,c(2,1,3:9)])
  ### and then run MigClim with "raster" type input.
  n<-MigClim.migrate (initDistrFile="InitialDist",
     hsMapFile="HSmap", barrierFile="Barrier",
     barrierType="weak", nrEnvChgSteps=5, nrDispSteps=5,
     dispKernel=c(1.0,0.4,0.16,0.06,0.03), initMatAge=1, fullMatAge=5,
     seedProdProb=c(0.01,0.08,0.5,0.92), lddFreq=0.1, minDist=6,
     maxDist=15, fullOutput=FALSE, simulName="mcTest", 
     overWrite=TRUE, testMode=FALSE, keepTempFiles=FALSE)}}
