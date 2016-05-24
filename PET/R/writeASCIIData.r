writeASCIIData <- function(data, filename, 
                  fileHeader=list(SignalDim=dim(data), 
                                  XYmin=c(-1,-1),
                                  DeltaXY=c(2/(nrow(data)-1), 2/(ncol(data)-1))),
                  Compress=TRUE, DebugLevel = "Normal"){

#########################################################
#
# The header of the file have the following form:
#
#     Description: char (optional)
#     SignalDim: int int
#     XYmin: double double
#     DeltaXY: double double
#
# After the header follow the raw data.
#
#########################################################
#
#
#  A <- matrix(1:6,nrow=2)
#  fileHeader <- list(SignalDim = c(2,3), XYmin = c(0,1), DeltaXY = c(1,1))
#  writeASCIIData(A ,"./test/NEUTEST2.dat", fileHeader)

 #######################################################
 #                  Checking the input-value
  DL1 <- logDebug(DebugLevel)
  DebugLevel <- DL1[[3]]
  DL2 <- DL1[[2]]
  DL1 <- DL1[[1]]
  
  if(!is.list(fileHeader))
      stop("'fileHeader' has to be of class 'list'.")
  lengthList <- length(fileHeader)
  if(lengthList != 3 && lengthList !=4)
      stop("The length of 'fileHeader' needs to be of length 3 or 4.")
  listNames <- names(fileHeader)
  if(lengthList == 3){
      if(!(all( listNames %in% c("SignalDim", "XYmin", "DeltaXY") ))){
          if(DL1) cat("WARNING: The names of 'fileHeader' doesn't correspond with defaults. \n")
          if(DL1) cat("It is assumed that the value in 'fileHeader' to be located in the correct order.\n")
          listNames <- FALSE } else listNames <- TRUE
  } else {
      if(!(all( listNames %in% c("Description", "SignalDim", "XYmin", 
                                 "DeltaXY") ))){
          if(DL1) cat("WARNING: The names of 'fileHeader' doesn't correspond with defaults. \n")
          if(DL1) cat("It is assumed that the value in 'fileHeader' to be located in the correct order.\n")
          listNames <- FALSE } else listNames <- TRUE
  }
  if(!(is.logical(Compress))){
      if(DL1) cat("WARNING: Parameter 'Compress' has to be TRUE or FALSE. Default is used. \n")
  }
 
 ########################################################
 #                  Main-Routine
  if(Compress){
      description <- paste(filename,".gz", sep="")
      if (DL1)
          cat("Open file: ", description, "\n")
      # Open connection
     con <- gzfile(description, "wb")
  } else {
      if (DL1)
          cat("Open file: ", filename, "\n")
      # Open connection
      con <- file(filename, "wb")
  }
  on.exit(close(con))
 
   
  if (DL2)
      cat("Creating the header of the file. \n")
  
  if(listNames){
     ###################################################
     #               Writing the header
     # writing 'Description' (optional), 'SignalDim', 
     # 'XYmin,' and 'DeltaXY'
      if(lengthList == 4){
          writeChar("Description", con, 11, eos=" ")
          write(fileHeader$Description, file=con)
      }
      writeChar("SignalDim", con, 9, eos=" ")
      write(fileHeader$SignalDim, file=con)
      writeChar("XYmin", con, 5, eos=" ")
      write(fileHeader$XYmin, file=con)
      writeChar("DeltaXY", con, 7, eos=" ")
      write(fileHeader$DeltaXY, file=con)
     ##################################################
     #             Writing the raw data
      if (DL2)
          cat("Writing the raw data.\n")
      index <- fileHeader$SignalDim
      dimIndex <- length(index)
      if (dimIndex > 2){
          perm.vec <- c(2,1, 3:dimIndex)
          write(x=aperm(data, perm.vec),
                file=con, ncolumns=index[2])
      } else if (dimIndex == 2){
          write(x=t(data), file=con, ncolumns=index[2])
      } else if (dimIndex == 1){
          write(x=data, file=con, ncolumns=index)
      } else
          stop("Dimension of the data is not supported.")
  } else {
     ###################################################
     #               Writing the header
     # writing 'Description' (optional), 'SignalDim', 
     # 'XYmin,' and 'DeltaXY'
      if(lengthList == 4){
          writeChar("Description", con, 11, eos=" ")
          write(fileHeader[[2]], file=con)
      }
      writeChar("SignalDim", con, 9, eos=" ")
      write(fileHeader[[lengthList-3]], file=con)
      writeChar("XYmin", con, 5, eos=" ")
      write(fileHeader[[lengthList-2]], file=con)
      writeChar("DeltaXY", con, 7, eos=" ")
      write(fileHeader[[lengthList-1]], file=con)
     ##################################################
     #             Writing the raw data
      if (DL2)
          cat("Writing the raw data.\n")
      index <- fileHeader[[lengthList-3]]
      dimIndex <- length(index)
      if (dimIndex > 2){
          perm.vec <- c(2,1, 3:dimIndex)
          write(x=aperm(data, perm.vec),
                file=con, ncolumns=index[2])
      } else if (dimIndex == 2){
          write(x=t(data), file=con, ncolumns=index[2])
      } else if (dimIndex == 1){
          write(x=data, file=con, ncolumns=index)
      } else
          stop("Dimension of the data is not supported.")
  }
  
  if (DL1)
      cat("Close file. \n")

  invisible(NULL)
}