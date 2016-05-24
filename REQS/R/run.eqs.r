run.eqs <- function(EQSpgm, EQSmodel, serial, Rmatrix = NA, datname = NA, LEN = 2000000)
{
  res <- call.eqs(EQSpgm = EQSpgm, EQSmodel = EQSmodel, serial = serial, Rmatrix = Rmatrix, datname = datname, LEN = LEN)
    
  if (!res) warning("EQS estimation not successful!")
  
  filedir.split <- strsplit(EQSmodel, "/")[[1]]
  n <- length(filedir.split)
  etsname <- strsplit(filedir.split[n], "\\.")[[1]][1]
  etsfile <- paste(etsname, ".ets",sep = "" )
  
  reslist <- read.eqs(etsfile)
  return(c(list(success = res),reslist))
}

