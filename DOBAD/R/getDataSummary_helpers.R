getDataSummary.CTMC_PO_many <- function(dat,file="dataSummary.rsav"){
  myIndepInts <- CTMCPO2indepIntervals.CTMC_PO_many(dat)
  myN <- length(dat@BDMCsPO);
  myIntLens <- myIndepInts[,3]
  myNumInts <- length(myIndepInts[,1]);
  myAvgIntLen <- mean(myIntLens);
  mySdIntLen <- sd(myIntLens)
  myMeanStartSt <- mean(myIndepInts[,1]);
  mySDstartSt <- sd(myIndepInts[,1]);
  myTotalT <- sum(myIndepInts[,3])
  ints.incr <- myIndepInts[myIndepInts[,2] > myIndepInts[,1],, drop=FALSE]
  ints.decr <- myIndepInts[myIndepInts[,2] < myIndepInts[,1],,drop=FALSE]
  ints.eq <- myIndepInts[myIndepInts[,2] == myIndepInts[,1],,drop=FALSE]
  myNumIntsI <- length(ints.incr[,1])
  myNumIntsD <- length(ints.decr[,1])
  myNumIntsE <- length(ints.eq[,1])
  myAvgI <-  mean(ints.incr[,2]-ints.incr[,1])
  myAvgD <- -mean(ints.decr[,2]-ints.decr[,1])
  ##myAvgE <- mean(ints.eq[,2]-ints.eq[,1]) ##0
  if (is.character(file) && file!=""){
    save(myIndepInts,
         myN,
         myIntLens,
         myNumInts,
         myAvgIntLen,
         mySdIntLen,
         myMeanStartSt,
         mySDstartSt, 
         myTotalT, 
         ints.incr,
         ints.decr,
         ints.eq,
         myNumIntsI,
         myNumIntsD,
         myNumIntsE,
         myAvgI, 
         myAvgD,
         file=file);
  }
  res <- list(myIndepInts=myIndepInts,
              myN=myN,
              myIntLens=myIntLens,
              myNumInts=myNumInts,
              myAvgIntLen=myAvgIntLen,
              mySdIntLen=mySdIntLen,
              myMeanStartSt=myMeanStartSt,
              mySDstartSt=mySDstartSt, 
              myTotalT=myTotalT, 
              ints.incr=ints.incr,
              ints.decr=ints.decr,
              ints.eq=ints.eq,
              myNumIntsI=myNumIntsI,
              myNumIntsD=myNumIntsD,
              myNumIntsE=myNumIntsE,
              myAvgI=myAvgI, 
              myAvgD=myAvgD)
  return(res)
}


## getDataSummary.CTMC_PO_many <- function(dat,file="dataSummary.rsav"){
##   myIndepInts <- CTMCPO2indepIntervals.CTMC_PO_many(dat)
##   myIntLens <- myIndepInts[,3]
##   myNumInts <- length(myIndepInts[,1]);
##   myAvgIntLen <- mean(myIntLens);
##   mySdIntLen <- sd(myIntLens)
##   myMeanStartSt <- mean(myIndepInts[,1]);
##   mySDstartSt <- sd(myIndepInts[,1]);
##   myTotalT <- sum(myIndepInts[,3])
##   ints.incr <- myIndepInts[myIndepInts[,2] > myIndepInts[,1],, drop=FALSE]
##   ints.decr <- myIndepInts[myIndepInts[,2] < myIndepInts[,1],,drop=FALSE]
##   ints.eq <- myIndepInts[myIndepInts[,2] == myIndepInts[,1],,drop=FALSE]
##   myNumIntsI <- length(ints.incr[,1])
##   myNumIntsD <- length(ints.decr[,1])
##   myNumIntsE <- length(ints.eq[,1])
##   myAvgI <-  mean(ints.incr[,2]-ints.incr[,1])
##   myAvgD <- -mean(ints.decr[,2]-ints.decr[,1])
##   ##myAvgE <- mean(ints.eq[,2]-ints.eq[,1])
##   save(myIndepInts,
##        myIntLens,
##        myNumInts,
##        myAvgIntLen,
##        mySdIntLen,
##        myMeanStartSt,
##        mySDstartSt, 
##        myTotalT, 
##        ints.incr,
##        ints.decr,
##        ints.eq,
##        myNumIntsI,
##        myNumIntsD,
##        myNumIntsE,
##        myAvgI, 
##        myAvgD,
##        file=file);  
## }
