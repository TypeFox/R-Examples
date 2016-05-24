call.eqs <- function(EQSpgm, EQSmodel, serial, Rmatrix = NA, datname = NA, LEN = 2000000)
{
  #EQSpgm ... EQS directory
  #EQSmodel ... path, where eqs-file is located
  #serial ... EQS serial number as character
  #LEN ... Number of working array unit.  By default, it is 2,000,000 8 bytes units
  
  #--------- change working directory where eqs file is located ---------
  if (!file.exists(EQSmodel)) stop("The .eqs file not found in the current folder!")
  
  filedir.split <- strsplit(EQSmodel, "/")[[1]]
  n <- length(filedir.split) 
  filedir <- paste(filedir.split[1:(n-1)], collapse = "/")
  if (n > 1) setwd(filedir)
   
  #--------- string specifications -----------
  outname <- strsplit(filedir.split[n], "\\.")[[1]][1]
  file.out <- paste(outname, ".out",sep = "" )
  
  lenstring <- paste("LEN=",as.integer(LEN), sep = "")
  filepathin <- paste("IN=", EQSmodel, sep = "")
  fileout <- paste("OUT=", file.out, sep = "")
  serstring <- paste("SER=", serial, "\n", sep = "")
  
  #----------- sanity check for input data ----------------
  if(length(Rmatrix) > 1) {                              #write data matrix, to file.dat, blank separated
    if(is.na(datname)) {
      warning(paste("No filename for data specified! ",outname,".dat is used", sep = "")) 
      datname <- paste(outname,".dat", sep = "")
    }
    write.table(as.matrix(Rmatrix), file = datname, col.names = FALSE, row.names = FALSE)
  }   
  
  
# IN=   input EQS command file (please note that input data is specified in the EQS command file)
# OUT=  EQS output file name

  EQScmd <- paste(deparse(EQSpgm), filepathin, fileout, lenstring, serstring)
  
  #RetCode <- system(EQScmd, intern = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL,
  #    show.output.on.console = FALSE, minimized = FALSE, invisible = FALSE)
      
  RetCode <- system(EQScmd, intern = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  if (RetCode == 0) {
    success <- TRUE
  } else {
    success <- FALSE
    warning("EQS estimation was not successful!") 
  }
  
  return(success = success)
}



