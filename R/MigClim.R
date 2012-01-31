#
# MigClim.R: The R functions for the MigClim package.
#
# Robin Engler & Wim Hordijk   Last modified: 25 January 2012
#


#
# MigClim.migrate: Initialize the MigClim method by writing the parameter
#                  values to file, and then run it.
#
MigClim.migrate <- function (initDistrFile="InitialDist", hsMapFile="HSmap",
                          barrierFile="", barrierType="weak", nrEnvChgSteps=1,
                          nrDispSteps=1, dispKernel=c(1.0,1.0), initMatAge=1,
                          fullMatAge=1, seedProdProb=c(1.0), rcThreshold=0,
                          lddFreq=0.0, minDist=NULL, maxDist=NULL,
                          fullOutput=FALSE, simulName="mySimul", 
                          overWrite=FALSE, testMode=FALSE, keepTempFiles=FALSE)
{
  #
  # If the user has entered a file name (as opposed to a dataframe or matrix)
  # then we emove any ".asc" or ".tif" extension that the user may have specified
  # in his/her filename.
  #
  if(is.character(initDistrFile))
  {
	  if (substr(initDistrFile, nchar(initDistrFile)-3, nchar(initDistrFile)) ==
	      ".asc")
	  {
	    initDistrFile <- strtrim(initDistrFile, nchar(initDistrFile)-4)
	  }
	  if (substr(hsMapFile, nchar(hsMapFile)-3, nchar(hsMapFile)) == ".asc")
	  {
	    hsMapFile <- strtrim(hsMapFile,nchar(hsMapFile)-4)
	  }
	  if (substr(initDistrFile, nchar(initDistrFile)-3, nchar(initDistrFile)) ==
	      ".tif")
	  {
	    initDistrFile <- strtrim(initDistrFile, nchar(initDistrFile)-4)
	  }
	  if (substr(hsMapFile, nchar(hsMapFile)-3, nchar(hsMapFile)) ==".tif")
	  {
	    hsMapFile <- strtrim(hsMapFile,nchar(hsMapFile)-4)
	  }
  }

  #
  # Detect the type of input given by the user. This can be any of
  # the following:
  #  -> datafraome or matrix
  #  -> ascii grid (.asc), geo-tiff (.tif)
  #     ESRI raster (no extension) or R raster (no extension).
  #
  if(require(raster, quietly=T)==F) stop("This function requires the 'raster' package. Please install 'raster' on your computer and try again.")
  RExt <- NA
  if(is.matrix(initDistrFile)) initDistrFile <- as.data.frame(initDistrFile)
  if(is.data.frame(initDistrFile)){
	  RExt <- ".DataFrame"
  } else{
	  if (file.exists(initDistrFile))
	  {
	    Rst <- try(raster(initDistrFile), silent=T)
	    if(class(Rst)[1]=="RasterLayer") RExt <- ""
	    rm(Rst)
	  }
	  if (file.exists(paste(initDistrFile,".tif",sep="")))
	  {
	    Rst <- try(raster(paste(initDistrFile,".tif",sep="")), silent=T)
	    if(class(Rst)[1]=="RasterLayer") RExt <- ".tif"
	    rm(Rst)
	    
	  }
	  if (file.exists(paste(initDistrFile,".asc",sep="")))
	  {
	    Rst <- try(raster(paste(initDistrFile,".asc",sep="")), silent=T)
	    if(class(Rst)[1]=="RasterLayer") RExt <- ".asc"
	    rm(Rst)
	  }
  }
  if (is.na(RExt))
  {
    stop ("Input data not recognized. Your input raster data must be in one of the following formats: ascii grid (.asc), geoTiff (.tif), ESRI grid or R raster (no extension). \n")
  }

  
  #
  # If the user chose to not allow overwriting of existing files
  # then we check that no such file already exists.
  if(overWrite==F){
	   
	  ### Check if output directory exists
	  if(file.exists(simulName)) stop("The output directory '", getwd(), "/", simulName, "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
	  
	  ### Check if any ".asc" files already exist.
	  if(RExt!=".asc"){
		  if(file.exists(paste(basename(initDistrFile),".asc",sep=""))) stop("The output directory '", getwd(), "/", paste(basename(initDistrFile),".asc",sep=""), "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
		  for(J in 1:nrEnvChgSteps) if(file.exists(paste(basename(hsMapFile), J,".asc",sep=""))) stop("The output directory '", getwd(), "/", paste(basename(hsMapFile), J,".asc",sep=""), "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
		  if (barrierFile!="") if(file.exists(paste(basename(barrierFile),".asc",sep=""))) stop("The output directory '", getwd(), "/", paste(basename(barrierFile),".asc",sep=""), "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
	  }
	  if(RExt==".asc"){
		  if(initDistrFile!=basename(initDistrFile)) if(file.exists(paste(basename(initDistrFile),".asc",sep=""))) stop("The output directory '", getwd(), "/", paste(basename(initDistrFile),".asc",sep=""), "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
		  if(hsMapFile!=basename(hsMapFile)) for(J in 1:nrEnvChgSteps) if(file.exists(paste(basename(hsMapFile), J,".asc",sep=""))) stop("The output directory '", getwd(), "/", paste(basename(hsMapFile), J,".asc",sep=""), "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
		  if(barrierFile!="") if(barrierFile!=basename(barrierFile)) if(file.exists(paste(basename(barrierFile),".asc",sep=""))) stop("The output directory '", getwd(), "/", paste(basename(barrierFile),".asc",sep=""), "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
	  }
	  if(RExt==".DataFrame"){
		  if(file.exists(paste(simulName, ".InitialDist.asc", sep=""))) stop("The output directory '", getwd(), "/", paste(simulName, ".InitialDist.asc", sep=""), "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
		  for(J in 1:nrEnvChgSteps) if(file.exists(paste(simulName, ".HSmap", J, ".asc", sep=""))) stop("The output directory '", getwd(), "/", paste(simulName, ".HSmap", J, ".asc", sep=""), "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
		  if (barrierFile!="") if(file.exists(paste(simulName, ".Barrier.asc", sep=""))) stop("The output directory '", getwd(), "/", paste(simulName, ".Barrier.asc", sep=""), "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")  
	  }
  }
  
  
  
  
  #
  # If the user has given the input as a matrix/dataframe, then we verify that
  # the data has the correct format. The correct format is as follows:
  # 'initDistrFile' = a data frame with 3 columns: X coordinate, Y coordinate and the species' initial distribution (0 or 1 values only).
  # 'hsMapFile' = a dataframe with ncol = nrEnvChgSteps. values must be in the range [0:1000]
  # 'barrierFile' = a dataframe or vector containing only values of 0 or 1.
  #
  if(RExt==".DataFrame")
  {
	cat("Converting data to ascii grid format... \n")  
	
	### Convert all input data to data frame objects.
	if(is.matrix(hsMapFile)) hsMapFile <- as.data.frame(hsMapFile)
	if(is.vector(hsMapFile)) hsMapFile <- as.data.frame(hsMapFile)
	if(is.character(barrierFile)) useBarrier <- FALSE else useBarrier <- TRUE
	if(is.matrix(barrierFile)) barrierFile <- as.data.frame(barrierFile)
	if(is.vector(barrierFile)) barrierFile <- as.data.frame(barrierFile)
	
	### Verify all inputs are of the same type
	if(!is.data.frame(hsMapFile)) stop("Data input error: the 'hsMapFile' data could not be converted to a dataframe. All inputs must be of the same type. \n")
	if(useBarrier) if(!is.data.frame(barrierFile)) stop("Data input error: the 'barrierFile' data could not be converted to a dataframe. all inputs must be of the same type. \n")
	
	### Verify all data frames have the correct number of rows and columns.
	if(ncol(initDistrFile)!=3) stop("Data input error. When entering 'initDistrFile' as a data frame or matrix, the data frame must have exactly 3 columns (in this order): X and Y coordinates, Initial distribution of the species. \n")
	if(ncol(hsMapFile)!=nrEnvChgSteps) stop("Data input error. When entering 'hsMapFile' as a data frame or matrix, the data frame must have a number of columns equal to nrEnvChgSteps. \n")
	if(nrow(hsMapFile)!=nrow(initDistrFile))  stop("Data input error. 'initDistrFile' and 'hsMapFile' must have the same number of rows.\n")
	if(useBarrier){
		if(ncol(barrierFile)!=1) stop("Data input error. When entering 'barrierFile' as a data frame, matrix or vector, the data must have a excatly 1 column. \n")
		if(nrow(barrierFile)!=nrow(initDistrFile))  stop("Data input error. 'initDistrFile' and 'barrierFile' must have the same number of rows.\n")
	}
	
	### Verify all data frames contain meaning full values.
	if(any(is.na(match(unique(initDistrFile[,3]), c(0,1))))) stop("Data input error: the 3rd column of 'initDistrFile' should contain only values of 0 or 1. \n")
	if(any(hsMapFile<0) | any(hsMapFile>1000)) stop("Data input error: all values in 'hsMapFile' must be in the range [0:1000]. \n")
	if(useBarrier){
		if(any(is.na(match(unique(barrierFile[,1]), c(0,1))))) stop("Data input error: 'barrierFile' should contain only values of 0 or 1. \n")
	}
	
	### Convert data frames to ascii grid files.
	if(require(SDMTools, quietly=T)==F) stop("This function requires the 'SDMTools' package. Please install 'SDMTools' on your computer and try again.")
	CreatedASCII <- paste(simulName, c("InitialDist.asc", paste("HSmap", 1:nrEnvChgSteps, ".asc", sep="")), sep=".")
	dataframe2asc(cbind(initDistrFile[,c(2,1,3)], hsMapFile), outdir=getwd(), filenames=CreatedASCII, gz=FALSE)
	if(useBarrier){
		dataframe2asc(cbind(initDistrFile[,c(2,1)], barrierFile), outdir=getwd(), filenames=paste(simulName, ".Barrier", sep=""), gz=FALSE)
		CreatedASCII <- c(CreatedASCII, paste(simulName, ".Barrier.asc", sep=""))
		barrierFile <- paste(simulName, ".Barrier", sep="")
	}
	initDistrFile <- paste(simulName, ".InitialDist", sep="")
	hsMapFile <- paste(simulName, ".HSmap", sep="")
	RExt <- ".asc"
  }
  
  
  #
  # Verify that all the input raster files do exist.
  #
  if (!file.exists(paste(initDistrFile,RExt,sep="")))
  {
    stop("'initDistrFile' could not be found.")
  }
  for (J in 1:nrEnvChgSteps)
  {
    if (!file.exists(paste(hsMapFile,J,RExt,sep="")))
    {
      stop("One or more of the 'hsMapFiles' could not be found.")
    }
  }
  if (barrierFile!="")
  {
    if (!file.exists(paste(barrierFile,RExt,sep="")))
    {
      stop("'barrierFile' could not be found.")
    }
  }
  
  #
  # Verify that parameters have meaningful values.
  #
  if (barrierFile!="")
  {
    if (!any(barrierType==c("weak","strong")))
    {
      stop("'barrierType' must be either 'weak' or 'strong'")
    }
  }
  if (!is.numeric(nrEnvChgSteps))
  {
    stop("'nrEnvChgSteps' must be a number in the range 1:99")
  }
  if (nrEnvChgSteps<1 | nrEnvChgSteps > 99)
  {
    stop("'nrEnvChgSteps' must be a number in the range 1:99")
  }
  if (!is.numeric(dispKernel))
  {
    stop("Values of 'dispKernel' must be numbers in the range ]0:1]")
  }
  if (any(dispKernel>1 | any(dispKernel<=0)))
  {
    stop("Values of 'dispKernel' must be numbers in the range ]0:1]")
  }
  if (!is.numeric(lddFreq))
  {
    stop("Data input error: 'lddFreq' must be a numeric value. \n")
  }
  if (lddFreq>0)
  {
    if (!is.numeric(minDist))
    {
      stop("Data input error: 'minDist' must be a numeric value. \n")
    }
    if (minDist <= length(dispKernel))
    {
      stop("Data input error: 'minDist' must be larger than the maximum distance given in the 'dispKernel'. \n")
    }
    if (!is.numeric(maxDist))
    {
      stop("Data input error: 'maxDist' must be a numeric value. \n")
    }
    if (maxDist < minDist)
    {
      stop("Data input error: 'maxDist' must be >= 'minDist'. \n")
    }
  } else
  {
    minDist <- maxDist <- 0
  }
  
  #
  # If the input format is not ascii grid, then we convert the files to
  # ascii grid format.
  # Note that we store the names of the created ascii files in the
  # "CreatedASCII" object.
  #
  if (RExt!=".asc")
  {
    cat("Converting data to ascii grid format... \n")
    Rst <- raster(paste(initDistrFile,RExt,sep=""))
    initDistrFile <- basename(initDistrFile)
    Rst2 <- writeRaster(Rst, filename=paste(initDistrFile,".asc",sep=""),
                        format="ascii", overwrite=TRUE, datatype="INT2S", NAflag=-9999)
    CreatedASCII <- paste(initDistrFile,".asc",sep="")
    for (J in 1:nrEnvChgSteps)
    {
      Rst <- raster(paste(hsMapFile,J,RExt,sep=""))
      Rst2 <- writeRaster(Rst, filename=paste(basename(hsMapFile),J,".asc",sep=""),
                          format="ascii", overwrite=TRUE, datatype="INT2S", NAflag=-9999)
      CreatedASCII <- c(paste(basename(hsMapFile),J,".asc",sep=""), CreatedASCII)
    }
    hsMapFile <- basename(hsMapFile)
    if (barrierFile!="")
    {
      Rst <- raster(paste(barrierFile,RExt,sep=""))
      barrierFile <- basename(barrierFile)
      Rst2 <- writeRaster(Rst, filename=paste(barrierFile,".asc",sep=""),
                          format="ascii", overwrite=TRUE, datatype="INT2S", NAflag=-9999)
      CreatedASCII <- c(paste(barrierFile,".asc",sep=""), CreatedASCII)
    }
    rm(Rst,Rst2)
  }
  
  #
  # If the input format is ascii grid, then we check that the files
  # are located in the working directory. If not, then we copy
  # the files to the working directory.
  if (RExt==".asc")
  {
    if(initDistrFile!=basename(initDistrFile)){
		file.copy(from=initDistrFile, to=basename(initDistrFile), overwrite=T)
		initDistrFile <- basename(initDistrFile)
		if(exists(CreatedASCII)) CreatedASCII <- c(paste(initDistrFile,".asc",sep=""), CreatedASCII) else CreatedASCII <- paste(initDistrFile,".asc",sep="")
    }
    if(hsMapFile!=basename(hsMapFile)){
		for(J in 1:nrEnvChgSteps){
			file.copy(from=paste(hsMapFile,J,sep=""), to=paste(basename(hsMapFile),J,sep=""), overwrite=T)
			if(exists(CreatedASCII)) CreatedASCII <- c(paste(basename(hsMapFile),J,".asc",sep=""), CreatedASCII) else CreatedASCII <- paste(basename(hsMapFile),J,".asc",sep="")
		}
		hsMapFile <- basename(hsMapFile)    
	}
    if(barrierFile!=""){
		if(barrierFile!=basename(barrierFile)){
			file.copy(from=barrierFile, to=basename(barrierFile), overwrite=T)
			barrierFile <- basename(barrierFile)
			if(exists(CreatedASCII)) CreatedASCII <- c(paste(barrierFile,".asc",sep=""), CreatedASCII) else CreatedASCII <- paste(barrierFile,".asc",sep="")
		}
    }
  }
  
  #
  # Verify that all raster have exactly the same dimensions and that they
  # contain apropriate values. "initDistrFile" and "barrierFile" should
  # contain only values of 0 or 1. "hsMapFile" should contain only values
  # in the range [0:1000]
  #
  Rst <- raster(paste(initDistrFile,".asc",sep=""))
  nrRows <- nrow(Rst)
  nrCols <- ncol(Rst)
  if(any(is.na(match(raster::unique(Rst), c(0,1)))))
  {
	  stop("Data input error: the 'initDistrFile' raster should contain only values of 0 or 1 \n")
  }
  for (J in 1:nrEnvChgSteps)
  {
    Rst <- raster(paste(hsMapFile,J,".asc",sep=""))
    if (nrow(Rst)!=nrRows | ncol(Rst)!=nrCols)
    {
      stop("Data input error: not all your rasters input data have the same dimensions. \n")
    }
    if (cellStats(Rst,"min")<0 | cellStats(Rst,"max")>1000)
    {
      stop("Data input error: not all habitat suitability rasters must have values in the range [0:1000]. \n")
    }
    rm(Rst)
  }
  if (barrierFile!="")
  {
    Rst <- raster(paste(barrierFile,".asc",sep=""))
    if (nrow(Rst)!=nrRows | ncol(Rst)!=nrCols)
    {
      stop("Data input error: not all your rasters input data have the same dimensions. \n")
    }
    if(any(is.na(match(raster::unique(Rst), c(0,1)))))
    {
      stop("Data input error: the 'barrierFile' raster should contain only values of 0 or 1 \n")
    }
    rm(Rst)
  }

  #
  # Create output directory.
  #
  if (file.exists(simulName)==T)
  {
    unlink(simulName, recursive=T)
  }
  if (dir.create(simulName)==F)
  {
    stop("unable to create a '", simulName,"'subdirectory in the current workspace. Make sure the '", simulName,"'subdirectory does not already exists and that you have write permission in the current workspace.")
  }
	
  #
  # Write the "simulName_params.txt" file to disk.
  #
  fileName <- paste(simulName, "/", simulName, "_params.txt", sep="")
  write (paste ("nrRows", nrRows), file=fileName, append=F)
  write (paste ("nrCols", nrCols), file=fileName, append=T)
  write (paste ("initDistrFile", initDistrFile), file=fileName, append=T)
  write (paste ("hsMapFile", hsMapFile), file=fileName, append=T)
  if (barrierFile != "")
  {
    write (paste ("barrierFile", barrierFile), file=fileName, append=T)
    write (paste ("barrierType", barrierType), file=fileName, append=T)
  }
  write (paste ("nrEnvChgSteps", nrEnvChgSteps), file=fileName, append=T)
  write (paste ("nrDispSteps", nrDispSteps), file=fileName, append=T)
  write (paste ("dispDist", length(dispKernel)), file=fileName, append=T)
  write (c ("dispKernel", dispKernel), file=fileName, append=T,
         ncolumns = length (dispKernel)+1)
  write (paste ("initMatAge", initMatAge), file=fileName, append=T)
  write (paste ("fullMatAge", fullMatAge), file=fileName, append=T)
  write (c ("seedProdProb", seedProdProb), file=fileName, append=T,
         ncolumns = length (seedProdProb)+1)
  write (paste ("rcThreshold", rcThreshold), file=fileName, append=T)
  if (lddFreq > 0.0)
  {
    write (paste ("lddFreq", lddFreq), file=fileName, append=T)
    write (paste ("minDist", minDist), file=fileName, append=T)
    write (paste ("maxDist", maxDist), file=fileName, append=T)
  }
  if (fullOutput){
    write ("fullOutput true", file=fileName, append=T)
  } else{
    write ("fullOutput false", file=fileName, append=T)
  }
  write (paste ("simulName", simulName), file=fileName, append=T)

  #
  # Call the C function.
  #
  if(!testMode){
	cat("Starting simulation for ", simulName, "...\n") 
    migrate <- .C("mcMigrate",
                  paste(simulName, "/", simulName, "_params.txt", sep=""),
                  nr=integer(1))
  }
  #
  # If ASCII grids were created in the MigClim.init() function, then we
  # delete them here (unless the user has set "keepTempFiles" to TRUE.
  #
  if(keepTempFiles) rm(CreatedASCII)
  if (exists("CreatedASCII"))
  {
    for (J in 1:length(CreatedASCII))
    {
      unlink(CreatedASCII[J])
    }
    rm(CreatedASCII)
  }
  
  #
  # If the user selected "testMode", then we delete the created ouput directory
  if(testMode) unlink(simulName, recursive=T)
  
  #
  # Return the number of output files created.
  #
  if(!testMode) return(migrate$nr)
  if(testMode){
    cat("Test for", simulName, "completed sucessfully.\n")  
    return(nrEnvChgSteps)
  } 
}





MigClim.plot <- function(asciiFile, outDir="", fileFormat="jpeg", fullOutput=FALSE){
	
	
	### Verify user input
	if(substr(asciiFile, nchar(asciiFile)-3, nchar(asciiFile))!=".asc") asciiFile <- paste(asciiFile, ".asc", sep="")
	if(fileFormat!="jpeg" & fileFormat!="png" & fileFormat!="inR") stop("Input error: 'fileFormat' must be one of 'jpeg' or 'png'.\n")
	if(!file.exists(asciiFile)) stop("Input error: ", asciiFile, " could not be found.\n")
	if(outDir!="") if(!file.exists(outDir)) stop("Input error: 'outDir' directory could not be found.\n")
	
	### load raster library
	if(require(raster, quietly=T)==F) stop("This function requires the 'raster' package. Please install 'raster' on your computer and try again.")
	
	### get root name and output directory
	rootName <- substr(basename(asciiFile), 1, nchar(basename(asciiFile))-11)
	if(outDir=="" & dirname(asciiFile)!=".") outDir <- dirname(asciiFile)
	outDir <- paste(outDir,"/",sep="")
	inDir <- paste(dirname(asciiFile),"/",sep="")
	if(inDir=="./") inDir <- ""
	
	### get all files
	if(fullOutput){
		stepName <- paste(rootName, "_step_", sep="")
		fileList <- list.files(dirname(asciiFile))
		fileList <- fileList[which(substr(fileList, 1, nchar(stepName))==stepName)]
		fileList <- c(basename(asciiFile), fileList)
		rm(stepName)
	} else fileList <- basename(asciiFile)
	fileList <- substr(fileList,1,nchar(fileList)-4) # remove ".asc" extension.
	
	for(fileName in fileList){
		
		#fileName <- fileList[1]
		cat("plotting data for", paste(fileName,".asc",sep=""), "\n")
		Rst <- raster(paste(inDir,fileName,".asc",sep=""))
		
		### Create color ramp
		rstVals <- sort(raster::unique(Rst))
		negativeNb <- length(which(rstVals<0))
		positiveNb <- length(which(rstVals>1 & rstVals<30000))
		zeroExists <- any(rstVals==0)
		oneExists <- any(rstVals==1)
		unilimtedExists <- any(rstVals==30000)
		Colors <- rep("yellow", negativeNb)
		if(zeroExists) Colors <- c(Colors, "white")
		if(oneExists) Colors <- c(Colors, "black")
		Colors <- c(Colors, rainbow(positiveNb, start=0, end=0.4))
		if(unilimtedExists) Colors <- c(Colors, "pink")
		
		### Create a graphic window that has the correct ratio (to match the map)
		#x11(width=7, height=7*((ymax(Rst)-ymin(Rst))/(xmax(Rst)-xmin(Rst))))
		#image(Rst, col=Colors, breaks=c(min(rstVals)-1,rstVals))
		#savePlot(filename= "mcTest_raster", type = "jpeg", device = dev.cur(), restoreConsole = TRUE)
		
		### Plot
		if(fileFormat=="jpeg"){
			jpeg(filename=paste(outDir,fileName,".jpg",sep=""), width = 2000, height = 2000*((ymax(Rst)-ymin(Rst))/(xmax(Rst)-xmin(Rst))), quality=90, res=300)
			image(Rst, col=Colors, breaks=c(min(rstVals)-1,rstVals))
			dev.off()
		}
		if(fileFormat=="png"){
			png(filename=paste(outDir,fileName,".png",sep=""), width = 2000, height = 2000*((ymax(Rst)-ymin(Rst))/(xmax(Rst)-xmin(Rst))), res=300)
			image(Rst, col=Colors, breaks=c(min(rstVals)-1,rstVals))
			dev.off()
		}
		if(fileFormat=="inR"){
			x11(width=7, height=7*((ymax(Rst)-ymin(Rst))/(xmax(Rst)-xmin(Rst))))
			image(Rst, col=Colors, breaks=c(min(rstVals)-1,rstVals))
		}
		
		### Release memory
		rm(Rst,rstVals,negativeNb,positiveNb,zeroExists,oneExists,unilimtedExists,Colors)
	}
	
}


#
# MigClim.genClust: Run the genetic clusters migration simulation.
#
MigClim.genClust <- function (hsMapFile="hsMap", barrierFile="barrier",
                              nrClusters=4, nrIterations=1, threshold=445,
                              outFile="out", initFile="")
{
  #
  # Get the number of rows and columns from the first input file.
  #
  Rst <- raster(paste(hsMapFile,"1.asc",sep=""))
  nrRows <- nrow(Rst)
  nrCols <- ncol(Rst)

  #
  # Call the genClust C function.
  #
  migrator <- .C("genClust", as.integer(nrRows), as.integer(nrCols),
                 as.integer(nrClusters), as.integer(nrIterations),
                 as.integer(threshold), hsMapFile, barrierFile, outFile,
                 initFile)
}


#
# MigClim.validate: Validate a genetic clusters migration output file.
#
MigClim.validate <- function (validateFile="Validation.txt", nrPoints=0,
                              simFile="out1.asc", nrClusters=4)
{
  #
  # Call the validation C function.
  #
  validate <- .C("validate", validateFile, as.integer(nrPoints), simFile,
                 as.integer(nrClusters), score=double(2))
  return (validate$score)
}
