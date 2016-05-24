#can be ctmcpo1 or ctmcpomany
CTMCPO2indepIntervals.CTMC_PO_1<- function(partialDat){
  numObs <- length(getStates(partialDat));
  startStates <- getStates(partialDat)[1:numObs-1];
  endStates <- getStates(partialDat)[2:numObs];
  timeInts <- getTimes(partialDat)[2:numObs]-getTimes(partialDat)[1:numObs-1]; #will all be the same often
  theArg <- matrix( c(startStates,endStates, timeInts), nrow=numObs-1,byrow=FALSE)    
}

CTMCPO2indepIntervals.CTMC_PO_many <- function(partialDat){
  res <- lapply(partialDat@BDMCsPO, CTMCPO2indepIntervals.CTMC_PO_1);
  n <- length(partialDat@BDMCsPO)
  startStates <- endStates<- integer();
  timeInts <- numeric();
  for (i in 1:n){
    startStates <- c(startStates, res[[i]][,1])
    endStates <- c(endStates, res[[i]][,2])
    timeInts <- c(timeInts, res[[i]][,3])
  }
  matrix(c(startStates,endStates,timeInts), ncol=3,byrow=FALSE);
}
