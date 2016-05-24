readSifData <- function(filename, DebugLevel="Normal")
{
# readSifData offer the possibility to read a sif-File that contain
# a sparse matrix.
# The sif-file has to be the following structure:
# At first: a Header ...
# At second: The raw data in the following form.
#     M integers n(m) that specify the numbers of value in each row m
#     Than follow alternating the nm indicies (integer) for the appropriate 
#     columns and nm accordingly values (float) for each row.
#
#
# source("readSifData.r")
# x <- readSifData("./XTest.sif.sif")
#

  DL1 <- logDebug(DebugLevel)
  DebugLevel <- DL1[[3]]
  DL2 <- DL1[[2]]
  DL1 <- DL1[[1]]

  if(!is.character(filename))
     stop("'filename' must be of type character.")
  if(length(filename) != 1)
     stop("Please specify exactly one 'filename'.")
  if(!file.exists(filename))
     stop("File '", filename, "' does not exist.")
  if(file.access(filename, 4))
     stop("No read permission for file ", filename)

 # Open connection
  con <- file(filename, "rb")
  on.exit(close(con)) # be careful ...  

 ##################################################################
 # Read the sparse-matrix
 #
  header <- list()
  header$XSamples <- readBin( con, integer(), n=1, size=4 )
  header$YSamples <- readBin( con, integer(), n=1, size=4 )
  header$DeltaX <- readBin( con, numeric(), n=1, size=4 )
  header$DeltaY <- readBin( con, numeric(), n=1, size=4 )
  header$Xmin <- readBin( con, numeric(), n=1, size=4 )
  header$Ymin <- readBin( con, numeric(), n=1, size=4 )
  header$ThetaSamples <- readBin( con, integer(), n=1, size=4 )
  header$RhoSamples <- readBin( con, integer(), n=1, size=4 )
  header$DeltaRho <- readBin( con, numeric(), n=1, size=4 )
  header$DeltaTheta <- readBin( con, numeric(), n=1, size=4 )
  header$ThetaMin <- readBin( con, numeric(), n=1, size=4 )
  header$RhoMin <- readBin( con, numeric(), n=1, size=4 )
  header$LowestALevel <- readBin( con, numeric(), n=1, size=4 )
  header$RadonKernel <- readBin( con, integer(), n=1, size=4 )
  header$OverSamp <- readBin( con, integer(), n=1, size=4 )
  header$Regularization <- readBin( con, numeric(), n=1, size=4 )
  header$Iterations <- readBin( con, integer(), n=1, size=4 )
  header$mode <- readBin( con, integer(), n=1, size=4 )
  header$UseFast <- readBin( con, integer(), n=1, size=4 )
  header$ConstrainMin <- readBin( con, numeric(), n=1, size=4 )
  header$ConstrainMax <- readBin( con, numeric(), n=1, size=4 )
  header$Alpha <- readBin( con, numeric(), n=1, size=4 )
  header$Beta <- readBin( con, numeric(), n=1, size=4 )
  header$IterationType <- readBin( con, integer(), n=1, size=4 )
  #header$SaveMatLab <- readBin( con, integer(), n=1, size=4 )

  header$KernelFileSave <- readBin( con, integer(), n=1, size=4 )
  header$SaveIterations <- readBin( con, integer(), n=1, size=4 )

  #StartFileName <- readBin( con, raw(), n=200 )
  #header$StartFileName <- formatC(rawToChar(StartFileName))

  RefFileName <- readBin( con, raw(), n=200 )
  header$RefFileName <- formatC(rawToChar(RefFileName))

  KernelFileName <- readBin( con, raw(), n=200 )
  header$KernelFileName <- formatC(rawToChar(KernelFileName))

  InFileName <- readBin( con, raw(), n=200 )
  header$InFileName <- formatC(rawToChar(InFileName))

  OutFileName <- readBin( con, raw(), n=200 )
  header$OutFileName <- formatC(rawToChar(OutFileName))

  SaveIterationsName <- readBin( con, raw(), n=200 )
  header$SaveIterationsName <- formatC(rawToChar(SaveIterationsName))

  header$M <- readBin( con, integer(), n=1, size=4 )
  header$N <- readBin( con, integer(), n=1, size=4 )

 ##################################################################
 # Read the sparse-matrix
 #
  MIndexNumber <- readBin( con, integer(), n=header$M, size=4 )
  SparseIndex <- list()
  SparseSignal <- list()

  pcType <- .Platform$OS.type # discern of '"unix"' or '"windows"'
  if (pcType=="unix") pcType <- 1
  else if (pcType=="windows") pcType <- 2
  else pcType <- 0

  if (DL1) 
    cat("Start of reading the sif-File of size", file.info(filename)$size,"Bytes \n")
  for(i in 1:header$M){
    if (pcType==1){
      if(DL1) cat("Progress:",trunc(i*100/header$M),"% \r");
      #Sys.sleep(1);   
    } else if (pcType==2) {
      if(DL1) cat("Progress:",trunc(i*100/header$M),"\r");
      flush.console();
      #Sys.sleep(1);
    }
    SparseIndex[[i]] <- readBin( con, integer(), n=MIndexNumber[i], size=4 )
    SparseSignal[[i]] <- readBin( con, numeric(), n=MIndexNumber[i], size=4 )
  }
  if  (pcType==1 || pcType==2) if(DL1) cat("\n")

 ##################################################################
 # Return reading values as list
 #
  invisible(list(MIndexNumber=MIndexNumber, 
                 SparseIndex=SparseIndex, 
                 SparseSignal=SparseSignal,
                 Header=header))

}