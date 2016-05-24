#
# The R functions for the MigClim package.
# Robin Engler & Wim Hordijk
# MigClim.migrate: Last modified: 11 May 2012
#
# MigClim.migrate: Initialize the MigClim method by writing the parameter
#                  values to file, and then run it.
#
MigClim.migrate <- function (iniDist="InitialDist", hsMap="HSmap", rcThreshold=0, 
                             envChgSteps=1, dispSteps=1, dispKernel=c(1.0,1.0), 
                             barrier="", barrierType="strong", 
                             iniMatAge=1, propaguleProd=c(1.0), 
                             lddFreq=0.0, lddMinDist=NULL, lddMaxDist=NULL,
                             simulName="MigClimTest", replicateNb=1, overWrite=FALSE,
                             testMode=FALSE, fullOutput=FALSE, keepTempFiles=FALSE)
{
  
  # Verify that the user has installed the "raster" and "SDMTools" library on his machine (this is no longer needed, R does this automatically).
  #if(require(raster, quietly=T)==F) stop("This function requires the 'raster' package. Please install 'raster' on your computer and try again.")
  #if(require(SDMTools, quietly=T)==F) stop("This function requires the 'SDMTools' package. Please install 'SDMTools' on your computer and try again.")

  # Verify that parameters have meaningful values.
  if(!is.numeric(rcThreshold)) stop("'rcThreshold' must be an integer number in the range [0:1000]. \n")
  if(rcThreshold<0 | rcThreshold > 1000) stop("'rcThreshold' must be an integer number in the range [0:1000]. \n")
  if(rcThreshold%%1!=0) stop("'rcThreshold' must be a number an integer number. \n")
  if(!is.numeric(envChgSteps)) stop("'envChgSteps' must be an integer number in the range [1:295]. \n")
  if(envChgSteps<1 | envChgSteps > 295) stop("'envChgSteps' must be an integer number in the range [1:295]. \n")
  if(envChgSteps%%1!=0) stop("'envChgSteps' must be a number an integer number. \n")
  if(!is.numeric(dispSteps)) stop("'dispSteps' must be a number in the range [1:99]. \n")
  if(dispSteps<1 | dispSteps > 99) stop("'dispSteps' must be a number in the range [1:99]. \n")
  if(dispSteps%%1!=0) stop("'dispSteps' must be a number an integer number. \n")
  
  if(!is.numeric(dispKernel)) stop("Values of 'dispKernel' must be numbers > 0 and <= 1. \n")
  if(any(dispKernel>1) | any(dispKernel<=0)) stop("Values of 'dispKernel' must be numbers > 0 and <= 1")
  if(barrier[1]!="") if(!any(barrierType==c("weak","strong"))) stop("'barrierType' must be either 'weak' or 'strong'. \n")
  
  if(!is.numeric(iniMatAge)) stop("'iniMatAge' must be an integer number > 0. \n")
  if(iniMatAge<=0 | iniMatAge%%1!=0) stop("'iniMatAge' must be an integer number > 0. \n")
  if(!is.numeric(propaguleProd)) stop("Values of 'propaguleProd' must be numbers > 0 and < 1. \n")
  if(any(propaguleProd>1) | any(propaguleProd<=0)) stop("Values of 'propaguleProd' must be numbers > 0 and <= 1. \n")
  if(length(propaguleProd)>1) if(propaguleProd[length(propaguleProd)]==1) stop("If the length of the 'propaguleProd' vector is > 1, then the last value cannot be 1. See the MigClim user guide 'MigClim.userGuide()' for detailed explanations on this paramter. \n")
  
  if(!is.numeric(lddFreq)) stop("Data input error: 'lddFreq' must be a numeric value. \n")
  if(lddFreq<0 | lddFreq>1) stop("'lddFreq' must be a number >= 0 and <= 1. \n")
  if(lddFreq>0){
    if(!is.numeric(lddMinDist)) stop("Data input error: 'lddMinDist' must be a numeric value. \n")
    if(!is.numeric(lddMaxDist)) stop("Data input error: 'lddMaxDist' must be a numeric value. \n")
    if(lddMinDist%%1!=0 | lddMaxDist%%1!=0) stop("'lddMinDist' and 'lddMaxDist' must be integer numbers. \n")
    if(lddMinDist <= length(dispKernel)) stop("Data input error: 'lddMinDist' must be larger than the length of the 'dispKernel'. \n")
    if(lddMaxDist < lddMinDist) stop("Data input error: 'lddMaxDist' must be >= 'lddMinDist'. \n")
  } else lddMinDist <- lddMaxDist <- 0
  
  if(!is.numeric(replicateNb)) stop("Data input error: 'replicateNb' must be a numeric, integer, value. \n")
  if(replicateNb<1 | replicateNb%%1!=0) stop("Data input error: 'replicateNb' must be an integer value >= 1. \n")
  
  if(!is.logical(overWrite)) stop("Data input error: 'overWrite' must be either TRUE or FALSE. \n")
  if(!is.logical(testMode)) stop("Data input error: 'testMode' must be either TRUE or FALSE. \n")
  if(!is.logical(fullOutput)) stop("Data input error: 'fullOutput' must be either TRUE or FALSE. \n")
  if(!is.logical(keepTempFiles)) stop("Data input error: 'keepTempFiles' must be either TRUE or FALSE. \n")
  
  if(!is.character(iniDist)) if(!is.matrix(iniDist) & !is.data.frame(iniDist)) stop("Data input error: 'iniDist' must be either a string, a data frame or a matrix. \n")
  if(!is.character(hsMap)) if(!is.matrix(hsMap) & !is.data.frame(hsMap) & !is.vector(hsMap)) stop("Data input error: 'hsMap' must be either a string, a data frame, a matrix or a vector. \n")
  if(!is.character(barrier)) if(!is.matrix(barrier) & !is.data.frame(barrier) & !is.vector(barrier)) stop("Data input error: 'barrier' must be either a string, a data frame, a matrix or a vector. \n")
  if(is.character(iniDist)){
	  # if "iniDist" is given in string format:
	  if(length(iniDist)>1 | length(hsMap)>1 | length(barrier)>1) stop("Data input error: When given as string input, 'iniDist', 'hsMap' and 'barrier' must have a length = 1.\n")
	  if(!is.character(hsMap)) stop("Data input error: 'iniDist' and 'hsMap' must have the same format: either both 'string' or both 'data frame/matrix/vector'. \n")
  	  if(barrier[1]!="") if(!is.character(barrier)) stop("Data input error: 'iniDist' and 'barrier' must have the same format: either both 'string' or both 'data frame/matrix/vector'. \n")
  }
  if(!is.character(iniDist)){
	  # if "iniDist" is given in data frame or matrix format:
	  if(is.character(hsMap)) stop("Data input error: 'iniDist' and 'hsMap' must have the same format: either both 'string' or both 'data frame/matrix/vector'. \n")
  	  if(barrier[1]!="") if(is.character(barrier)) stop("Data input error: 'iniDist' and 'barrier' must have the same format: either both 'string' or both 'data frame/matrix/vector'. \n") 
  }
  
  
  
  # If the user has entered a file name (as opposed to a dataframe or matrix) then we remove
  # any ".asc" or ".tif" extension that the user may have specified in his/her filename.
  #
  if(is.character(iniDist)){
	  if (substr(iniDist, nchar(iniDist)-3, nchar(iniDist)) == ".asc") iniDist <- strtrim(iniDist, nchar(iniDist)-4)
	  if (substr(hsMap, nchar(hsMap)-3, nchar(hsMap)) == ".asc") hsMap <- strtrim(hsMap,nchar(hsMap)-4)
	  if (substr(iniDist, nchar(iniDist)-3, nchar(iniDist)) == ".tif") iniDist <- strtrim(iniDist, nchar(iniDist)-4)
	  if (substr(hsMap, nchar(hsMap)-3, nchar(hsMap)) ==".tif") hsMap <- strtrim(hsMap,nchar(hsMap)-4)
  }

  # Detect the type of input given by the user. This can be any of the following:
  #  -> dataframe or matrix
  #  -> ascii grid (.asc), geo-tiff (.tif),
  #     ESRI raster (no extension), or R raster (no extension).
  #
  RExt <- NA
  if(is.matrix(iniDist)) iniDist <- as.data.frame(iniDist)  #if the user input is a matrix, we convert it to a data frame.
  if(is.data.frame(iniDist)){
	  RExt <- ".DataFrame"
  } else{
	  if (file.exists(iniDist))
	  {
	    Rst <- try(raster(iniDist), silent=T)
	    if(class(Rst)[1]=="RasterLayer") RExt <- ""
	    rm(Rst)
	  }
	  if (file.exists(paste(iniDist,".tif",sep="")))
	  {
	    Rst <- try(raster(paste(iniDist,".tif",sep="")), silent=T)
	    if(class(Rst)[1]=="RasterLayer") RExt <- ".tif"
	    rm(Rst)
	    
	  }
	  if (file.exists(paste(iniDist,".asc",sep="")))
	  {
	    Rst <- try(raster(paste(iniDist,".asc",sep="")), silent=T)
	    if(class(Rst)[1]=="RasterLayer") RExt <- ".asc"
	    rm(Rst)
	  }
  }
  if(is.na(RExt)) stop ("Input data not recognized. Your input raster data must be in one of the following formats: ascii grid (.asc), geoTiff (.tif), ESRI grid or R raster (no extension). \n")


  
  # If the user chose to not allow overwriting of existing files (overWrite==F)
  # then we check that no future output file already exists.
  if(overWrite==F){
	   
	  ### Check if output directory exists
	  if(file.exists(simulName)) stop("The output directory '", getwd(), "/", simulName, "' already exists. \n Delete this directory or set 'overWrite=TRUE' in the function's parameters.\n")
	  
	  ### Check if any ".asc" files already exist.
	  if(RExt!=".asc"){
		  if(file.exists(paste(basename(iniDist),".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(iniDist),".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
		  for(J in 1:envChgSteps) if(file.exists(paste(basename(hsMap), J,".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(hsMap), J,".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
		  if (barrier!="") if(file.exists(paste(basename(barrier),".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(barrier),".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
	  }
	  if(RExt==".asc"){
		  if(iniDist!=basename(iniDist)) if(file.exists(paste(basename(iniDist),".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(iniDist),".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
		  if(hsMap!=basename(hsMap)) for(J in 1:envChgSteps) if(file.exists(paste(basename(hsMap), J,".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(hsMap), J,".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
		  if(barrier!="") if(barrier!=basename(barrier)) if(file.exists(paste(basename(barrier),".asc",sep=""))) stop("The output file '", getwd(), "/", paste(basename(barrier),".asc",sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
	  }
	  if(RExt==".DataFrame"){
		  if(file.exists(paste(simulName, ".InitialDist.asc", sep=""))) stop("The output file '", getwd(), "/", paste(simulName, ".InitialDist.asc", sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
		  for(J in 1:envChgSteps) if(file.exists(paste(simulName, ".HSmap", J, ".asc", sep=""))) stop("The output file '", getwd(), "/", paste(simulName, ".HSmap", J, ".asc", sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")
		  if (barrier!="") if(file.exists(paste(simulName, ".Barrier.asc", sep=""))) stop("The output file '", getwd(), "/", paste(simulName, ".Barrier.asc", sep=""), "' already exists. \n Delete this file or set 'overWrite=TRUE' in the function's parameters.\n")  
	  }
  }
  
  
   
  # If the user has given the input as a matrix/dataframe, then we verify that
  # the data has the correct format. The correct format is as follows:
  # 'iniDist' = a data frame with 3 columns: X coordinate, Y coordinate and the species' initial distribution (0 or 1 values only).
  # 'hsMap' = a dataframe with ncol = envChgSteps. values must be in the range [0:1000]
  # 'barrier' = a dataframe or vector containing only values of 0 or 1.
  #
  if(RExt==".DataFrame")
  {
	cat("Converting data to ascii grid format... \n")  
	
	### Convert all input data to data frame objects.
	if(is.matrix(hsMap)) hsMap <- as.data.frame(hsMap)
	if(is.vector(hsMap)) hsMap <- as.data.frame(hsMap)
	if(is.character(barrier)) useBarrier <- FALSE else useBarrier <- TRUE
	if(is.matrix(barrier)) barrier <- as.data.frame(barrier)
	if(is.vector(barrier)) barrier <- as.data.frame(barrier)
	
	### Verify all inputs are of the data frame type (note: iniDist has already been checked earlier).
	if(!is.data.frame(hsMap)) stop("Data input error: the 'hsMap' data could not be converted to a dataframe. All inputs must be of the same type. \n")
	if(useBarrier) if(!is.data.frame(barrier)) stop("Data input error: the 'barrier' data could not be converted to a dataframe. all inputs must be of the same type. \n")
	
	### Verify all data frames have the correct number of rows and columns.
	if(ncol(iniDist)!=3) stop("Data input error. When entering 'iniDist' as a data frame or matrix, the data frame must have exactly 3 columns (in this order): X and Y coordinates, Initial distribution of the species. \n")
	if(ncol(hsMap)!=envChgSteps) stop("Data input error. When entering 'hsMap' as a data frame or matrix, the data frame must have a number of columns equal to envChgSteps. \n")
	if(nrow(hsMap)!=nrow(iniDist))  stop("Data input error. 'iniDist' and 'hsMap' must have the same number of rows.\n")
	if(useBarrier){
		if(ncol(barrier)!=1) stop("Data input error. When entering 'barrier' as a data frame, matrix or vector, the data must have a excatly 1 column. \n")
		if(nrow(barrier)!=nrow(iniDist))  stop("Data input error. 'iniDist' and 'barrier' must have the same number of rows.\n")
	}
	
	### Verify all data frames contain meaningful values.
	if(any(is.na(match(unique(iniDist[,3]), c(0,1))))) stop("Data input error: the 3rd column of 'iniDist' should contain only values of 0 or 1. \n")
	if(any(hsMap<0) | any(hsMap>1000)) stop("Data input error: all values in 'hsMap' must be in the range [0:1000]. \n")
	if(useBarrier) if(any(is.na(match(unique(barrier[,1]), c(0,1))))) stop("Data input error: 'barrier' should contain only values of 0 or 1. \n")
	
	### Convert data frames to ascii grid files.
	CreatedASCII <- paste(simulName, c("InitialDist.asc", paste("HSmap", 1:envChgSteps, ".asc", sep="")), sep=".")
	dataframe2asc(cbind(iniDist[,c(2,1,3)], hsMap), outdir=getwd(), filenames=CreatedASCII, gz=FALSE)
	if(useBarrier){
		dataframe2asc(cbind(iniDist[,c(2,1)], barrier), outdir=getwd(), filenames=paste(simulName, ".Barrier", sep=""), gz=FALSE)
		CreatedASCII <- c(CreatedASCII, paste(simulName, ".Barrier.asc", sep=""))
		barrier <- paste(simulName, ".Barrier", sep="")
	}
	iniDist <- paste(simulName, ".InitialDist", sep="")
	hsMap <- paste(simulName, ".HSmap", sep="")
	RExt <- ".asc"
  }
  
  
  # Verify that all the input raster files do exist.
  if(!file.exists(paste(iniDist,RExt,sep=""))) stop(paste("The 'iniDist' file '", iniDist, RExt, "' could not be found.\n", sep=""))
  for(J in 1:envChgSteps){
    if(!file.exists(paste(hsMap,J,RExt,sep=""))) stop(paste("The 'hsMap' file '", hsMap, J, RExt, "' could not be found.\n",
                                                            "The naming convention for hsMap files is 'hsMap basename + 1', 'hsMap basename + 2', etc...\n",
                                                            "e.g. if you set 'hsMap='habitatSuitMap'' then your first hsMap file must be named 'habitatSuitMap1'.\n",
                                                            "the following hsMap file must be named 'habitatSuitMap2', 'habitatSuitMap3' and so on.\n", sep=""))
  }
  if(barrier!="") if(!file.exists(paste(barrier,RExt,sep=""))) stop(paste("The 'barrier' file '", barrier, RExt, "' could not be found.\n", sep=""))

  
  # If the input format is not ascii grid, then we convert the files to ascii grid format.
  # Note that we store the names of the created ascii files in the "CreatedASCII" object.
  #
  if (RExt!=".asc"){
    cat("Converting data to ascii grid format... \n")
    Rst <- raster(paste(iniDist,RExt,sep=""))
    iniDist <- basename(iniDist)
    Rst2 <- writeRaster(Rst, filename=paste(iniDist,".asc",sep=""), format="ascii", overwrite=TRUE, datatype="INT2S", NAflag=-9999)
    CreatedASCII <- paste(iniDist,".asc",sep="")
    for (J in 1:envChgSteps){
      Rst <- raster(paste(hsMap,J,RExt,sep=""))
      Rst2 <- writeRaster(Rst, filename=paste(basename(hsMap),J,".asc",sep=""), format="ascii", overwrite=TRUE, datatype="INT2S", NAflag=-9999)
      CreatedASCII <- c(paste(basename(hsMap),J,".asc",sep=""), CreatedASCII)
    }
    hsMap <- basename(hsMap)
    if (barrier!=""){
      Rst <- raster(paste(barrier,RExt,sep=""))
      barrier <- basename(barrier)
      Rst2 <- writeRaster(Rst, filename=paste(barrier,".asc",sep=""), format="ascii", overwrite=TRUE, datatype="INT2S", NAflag=-9999)
      CreatedASCII <- c(paste(barrier,".asc",sep=""), CreatedASCII)
    }
    rm(Rst,Rst2)
  }
  
  
  # If the input format is ascii grid, then we check that the files are located 
  # in the working directory. If not, we copy the files to the working directory.
  if (RExt==".asc"){
    if(iniDist!=basename(iniDist)){
		file.copy(from=paste(iniDist,".asc",sep=""), to=paste(basename(iniDist),".asc",sep=""), overwrite=T)
		iniDist <- basename(iniDist)
		if(exists("CreatedASCII")) CreatedASCII <- c(paste(iniDist,".asc",sep=""), CreatedASCII) else CreatedASCII <- paste(iniDist,".asc",sep="")
    }
    if(hsMap!=basename(hsMap)){
		for(J in 1:envChgSteps){
			file.copy(from=paste(hsMap,J,".asc",sep=""), to=paste(basename(hsMap),J,".asc",sep=""), overwrite=T)
			if(exists("CreatedASCII")) CreatedASCII <- c(paste(basename(hsMap),J,".asc",sep=""), CreatedASCII) else CreatedASCII <- paste(basename(hsMap),J,".asc",sep="")
		}
		hsMap <- basename(hsMap)    
	}
    if(barrier!=""){
		if(barrier!=basename(barrier)){
			file.copy(from=paste(barrier,".asc",sep=""), to=paste(basename(barrier),".asc",sep=""), overwrite=T)
			barrier <- basename(barrier)
			if(exists("CreatedASCII")) CreatedASCII <- c(paste(barrier,".asc",sep=""), CreatedASCII) else CreatedASCII <- paste(barrier,".asc",sep="")
		}
    }
  }
  
  
  # Verify that all ascii grid files have a correct structure and that their NoData value (if any) is set to -9999
  #
  noDataVal <- getNoDataValue(paste(iniDist,".asc",sep=""))
  if(!is.na(noDataVal)){
    if(noDataVal=="ErrorInFile") stop("Data input error: the 'iniDist' ascii grid file does not have the correct structure.\n")
    if(noDataVal!=-9999) stop("Data input error: the 'iniDist' ascii grid file must have 'NoData' values set to -9999.\n")
  }
  for(J in 1:envChgSteps){
	noDataVal <- getNoDataValue(paste(hsMap,J,".asc",sep=""))
    if(!is.na(noDataVal)){
      if(noDataVal=="ErrorInFile") stop("Data input error: one or more 'hsMap' ascii grid files do not have the correct structure.\n")
      if(noDataVal!=-9999) stop("Data input error: all 'hsMap' ascii grid files must have 'NoData' values set to -9999.\n")
    }    
  }
  if(barrier!=""){
    noDataVal <- getNoDataValue(paste(barrier,".asc",sep=""))
    if(!is.na(noDataVal)){
      if(noDataVal=="ErrorInFile") stop("Data input error: the 'Barrier' ascii grid file does not have the correct structure.\n")
      if(noDataVal!=-9999) stop("Data input error: the 'Barrier' ascii grid file must have 'NoData' values set to -9999.\n")
    }
  }
  rm(noDataVal)
  
  
  # Verify that all raster have exactly the same dimensions and that they contain apropriate values. 
  # "iniDist" and "barrier" should contain only values of 0 or 1. "hsMap" should contain only values in the range [0:1000].
  #
  Rst <- raster(paste(iniDist,".asc",sep=""))
  nrRows <- nrow(Rst)
  nrCols <- ncol(Rst)
  if(any(is.na(match(raster::unique(Rst), c(0,1))))) stop("Data input error: the 'iniDist' raster should contain only values of 0 or 1. \n")
  #if(dataType(Rst)!="INT2U" & dataType(Rst)!="INT1U") stop("Data input error: the 'iniDist' layer must contain integer values (8 or 16-bit unsigned integers). The R 'dataType' code for 8-bit and 16-bit unsigned integers is 'INT1U' and 'INT2U'.")
  for(J in 1:envChgSteps){
    Rst <- raster(paste(hsMap,J,".asc",sep=""))
    #if(dataType(Rst)!="INT2U") stop("Data input error: all habitat suitability rasters must contain integer values (16-bit unsigned integers) in the range 0 to 1000. The R 'dataType' code for 16-bit unsigned integers is 'INT2U'.")
    if(nrow(Rst)!=nrRows | ncol(Rst)!=nrCols) stop("Data input error: not all your rasters input data have the same dimensions. \n")
    if(cellStats(Rst,"min")<0 | cellStats(Rst,"max")>1000) stop("Data input error: all habitat suitability rasters must have values in the range [0:1000]. \n")
    rm(Rst)
  }
  if(barrier!=""){
    Rst <- raster(paste(barrier,".asc",sep=""))
    #if(dataType(Rst)!="INT2U" & dataType(Rst)!="INT1U") stop("Data input error: the 'barrier' layer must contain integer values (8 or 16-bit unsigned integers). The R 'dataType' code for 8-bit and 16-bit unsigned integers is 'INT1U' and 'INT2U'.")
    if(nrow(Rst)!=nrRows | ncol(Rst)!=nrCols) stop("Data input error: not all your rasters input data have the same dimensions. \n")
    if(any(is.na(match(raster::unique(Rst), c(0,1))))) stop("Data input error: the 'barrier' raster should contain only values of 0 or 1. \n")
    rm(Rst)
  }

    
	  
  # Create output directory.
  if (file.exists(simulName)==T) unlink(simulName, recursive=T)
  if (dir.create(simulName)==F) stop("unable to create a '", simulName,"'subdirectory in the current workspace. Make sure the '", simulName,"'subdirectory does not already exists and that you have write permission in the current workspace.\n")
  
  # Write the "simulName_params.txt" file to disk.
  fileName <- paste(simulName, "/", simulName, "_params.txt", sep="")
  write(paste("nrRows", nrRows), file=fileName, append=F)
  write(paste("nrCols", nrCols), file=fileName, append=T)
  write(paste("iniDist", iniDist), file=fileName, append=T)
  write(paste("hsMap", hsMap), file=fileName, append=T)
  write(paste("rcThreshold", rcThreshold), file=fileName, append=T)
  write(paste("envChgSteps", envChgSteps), file=fileName, append=T)
  write(paste("dispSteps", dispSteps), file=fileName, append=T)
  write(paste("dispDist", length(dispKernel)), file=fileName, append=T)
  write(c("dispKernel", dispKernel), file=fileName, append=T, ncolumns=length(dispKernel)+1)
  if(barrier!=""){
    write(paste("barrier", barrier), file=fileName, append=T)
    write(paste("barrierType", barrierType), file=fileName, append=T)
  }
  write(paste("iniMatAge", iniMatAge), file=fileName, append=T)
  write(paste("fullMatAge", iniMatAge + length(propaguleProd)), file=fileName, append=T)
  write(c("propaguleProd", propaguleProd), file=fileName, append=T, ncolumns=length(propaguleProd)+1)
  if(lddFreq > 0.0){
    write(paste("lddFreq", lddFreq), file=fileName, append=T)
    write(paste("lddMinDist", lddMinDist), file=fileName, append=T)
    write(paste("lddMaxDist", lddMaxDist), file=fileName, append=T)
  }
  if(fullOutput) write("fullOutput true", file=fileName, append=T) else write("fullOutput false", file=fileName, append=T)
  write(paste("replicateNb", replicateNb), file=fileName, append=T)
  write(paste("simulName", simulName), file=fileName, append=T)
  
  
  # Call the C function.
  if(!testMode){
	cat("Starting simulation for ", simulName, "...\n") 
    migrate <- .C("mcMigrate", paste(simulName, "/", simulName, "_params.txt", sep=""), nr=integer(1))
  }

	  
  # If ASCII grids were created in the MigClim.init() function, then we
  # delete them here (unless the user has set "keepTempFiles = TRUE".
  if(keepTempFiles) rm(CreatedASCII)
  if(exists("CreatedASCII")){
    for (J in 1:length(CreatedASCII)) unlink(CreatedASCII[J])
    rm(CreatedASCII)
  }
  
  
  # If the user has set replicateNb > 1 then we generate a final, averaged, output.
  # The individual outputs are conserved, though.
  if(replicateNb>1 & !testMode){
	  
	  #Average the "_stats.txt" files 
	  statsFile <- read.table(paste(simulName,"/",simulName,"1_stats.txt",sep=""), header=T, as.is=T)
	  for(J in 2:replicateNb){
		  tempFile <- read.table(paste(simulName,"/",simulName, J,"_stats.txt",sep=""), header=T, as.is=T)
		  statsFile <- statsFile + tempFile
		  rm(tempFile)
	  }
	  statsFile <- round(statsFile/replicateNb, 2)
	  write.table(statsFile, file=paste(simulName,"/",simulName,"_stats.txt",sep=""), quote=F, row.names=F, sep="\t")
	  
	  #Average the "_summary.txt" files 
	  statsFile <- read.table(paste(simulName,"/",simulName,"1_summary.txt",sep=""), header=T, as.is=T)
	  for(J in 2:replicateNb){
		  tempFile <- read.table(paste(simulName,"/",simulName, J,"_summary.txt",sep=""), header=T, as.is=T)
		  statsFile <- rbind(statsFile, tempFile)
		  rm(tempFile)
	  }
	  statsFile[replicateNb+1,1] <- simulName
	  statsFile[replicateNb+1,2:ncol(statsFile)] <- round(apply(statsFile[1:replicateNb,2:ncol(statsFile)], 2, mean), 2)
	  write.table(statsFile, file=paste(simulName,"/",simulName,"_summary.txt",sep=""), quote=F, row.names=F, sep="\t")
	  rm(statsFile)
  }

    
  # If the user selected "testMode", then we delete the created ouput directory
  if(testMode) unlink(simulName, recursive=T)
  
  # Return the number of output files created.
  if(!testMode){
	if(migrate$nr==envChgSteps) cat("Simulation ", simulName, " completed successfully. Outputs stored in ", getwd(), "/", simulName,"\n", sep="")  
    return(migrate$nr)
  }
  if(testMode){
    cat("Test for", simulName, "completed sucessfully.\n")  
    return(envChgSteps)
  } 
}




###
### This function checks the structure of an ESRI ascii grid file and, if the structure is correct,
### returns its NoData value.
### If the structure of the file is not correct, the function returns a string: "ErrorInFile"
### If no NoData value is indicated (which is possible, since this information is optional),
### the function returns NA.
### 
getNoDataValue <- function(fileName){
	
	# verify that the ascii file has the correct structure:	
	noDataVal <- "ErrorInFile"
	fileStruct <- c("ncols","nrows","xllcorner","yllcorner","cellsize")
	fileStruct2 <- c(5,5,9,9,8)
	for(J in 0:4){
		lineVal <- scan(file=fileName, what="character", nlines=1, skip=J, quiet=TRUE)
		if(length(lineVal)!=2) return(noDataVal)
		if(nchar(lineVal[1])!=fileStruct2[J+1]) return(noDataVal)
		if(length(grep(fileStruct[J+1], lineVal[1], ignore.case=T))!=1) return(noDataVal)
	}
	
	# get NoData value
	# (note that this line is optional in the file. If the line is missing we return NA.
	lineVal <- scan(file=fileName, what="character", nlines=1, skip=5, quiet=TRUE)
	if(length(lineVal)<2) return(noDataVal)
	noDataVal <- NA
	if(length(grep("NODATA_value", lineVal[1], ignore.case=T))==1 & nchar(lineVal[1])==12) noDataVal <- lineVal[2]

	return(noDataVal)
}


