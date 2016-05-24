readASCIIData <- function(filename, DebugLevel = "Normal"){

#########################################################
#
# The header of the file have to be the following form:
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
#  B <- readASCIIData("./test/NEUTEST.dat")
#

  DL1 <- logDebug(DebugLevel)
  DebugLevel <- DL1[[3]]
  DL2 <- DL1[[2]]
  DL1 <- DL1[[1]]
  
  sep <- ""
  
  if (DL1)
      cat("Open file: ",filename, "\n")
  
  # extraction of file-extension
  fiEx <- unlist(strsplit(filename, ".", fixed=TRUE))
  fiEx <- fiEx[length(fiEx)]

  # Open connection
  if (fiEx=="gz")
      con <- gzfile(filename, "rb")
  else
      con <- file(description=filename, "rb")
  
  on.exit(close(con)) # be careful ...
  
 #######################################################
 #                  Reading header
 # reading 'Description' and 'SignalDim'  
  if (DL2)
      cat("Reading the header of the file.\n")
  temp <- try(scan(con, what = "", n=1, sep=sep, quiet=TRUE))
  if (class(temp) == "try-error")
       stop("File read error. Perhaps the file doesn't come up to the necessary form.")
  if (substr(temp,1,10) == "SignalDim"){
      ans <- list(SignalDim = NULL, XYmin = NULL, DeltaXY = NULL)
      ans$SignalDim <- scan(con, nlines=1, sep=sep, quiet=TRUE)
  } else if (substr(temp,1,12) == "Description"){
      ans <- list(Description = NULL, SignalDim = NULL, XYmin = NULL, 
                  DeltaXY = NULL)
      ans$Description <- readLines(con, n=1)
      tmp <- scan(con, what = "", n=1, sep=sep, quiet=TRUE)
      if (tmp != "SignalDim")
          stop("Unexcept file header. Parameter ",tmp," in second line not allowed.")
      ans$SignalDim <- scan(con, nlines=1, sep=sep, quiet=TRUE)
  } else
      stop("Unexcept file header. First line doesn't correspond with the file
            format.")

  # reading	'XYmin'
  tmp <- scan(con, what = "", n=1, sep=sep, quiet=TRUE)
      if (tmp != "XYmin")
          stop("Unexcept file header. Parameter ",tmp," doesn't correspond with 
                the file format.")
  ans$XYmin <- scan(con, nlines=1, sep=sep, quiet=TRUE)

  # reading	'DeltaXY'
  tmp <- scan(con, what = "", n=1, sep=sep, quiet=TRUE)
      if (tmp != "DeltaXY")
          stop("Unexcept file header. Parameter ",tmp," doesn't correspond with 
                the file format..")
  ans$DeltaXY <- scan(con, nlines=1, sep=sep, quiet=TRUE)

 
 #######################################################
 #                  Reading the raw data
  if (DL2)
      cat("Reading the raw data.\n")
  index <- ans$SignalDim
  dimIndex <- length(index)
  if (dimIndex > 2){
      perm.vec <- c(2,1, 3:dimIndex)
      Signal <- aperm(array(scan(con, sep=sep, quiet=TRUE),
                          dim=index[perm.vec]),perm.vec)
  } else if (dimIndex == 2){
      Signal <- matrix(scan(con, sep=sep, quiet=TRUE), nrow=index[1], 
                       ncol=index[2], byrow=TRUE)
  } else if (dimIndex == 1){
      Signal <- matrix(scan(con, sep=sep, quiet=TRUE), nrow=1, ncol=index, 
                       byrow=TRUE)
  } else
      stop("The dimension of data is not supported.")
  
  if (DL1)
      cat("Close file. \n")
  
  return(list(Signal=Signal, Header=ans ))
}
