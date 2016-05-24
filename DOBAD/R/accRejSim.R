##Code for accept-reject simulation.  This code is not optimal and should not be
## exported.


############Code for sampling marginally and picking out conditionally


#accept-reject algorithm
#Count for only those sims that match data. Returns vec of T/F w same length as sims
condIndices <- function(sims,  data){
  if (is(data, "CTMC_PO_1")){
    dataTimes <- getTimes(data);
    dataStates <- getStates(data);
  }
  else {
    dataTimes <- data$times;
    dataStates <- data$states;
  }
  waits <- NULL;  jumps <- NULL;
  acceptables <- 0;
  Holdtime <- Nplus <- Nminus <- 0;
  simsStates <- t(sapply(sims, function(x){getMCstate(CTMC=x, Ts=dataTimes)}));
  indices <- apply(simsStates, 1, function(x){identical(x, dataStates)});
  indices; 
}


###Returns indices of the sims that have the passed-in summary statistics
###Sims is a list of componenets of class "CTMC"
### epsilon is the tolerance for holdtimes.
### MatchingChains it he complete data to be conditioned on
indicesCondSuffStats <- function(sims, matchingChain, eps){
  ARfn <- function(mySim){ #"accept-reject" function
    dataCounts <- NijBD(matchingChain);
    dataHoldtime <- BDsummaryStats(matchingChain)["Holdtime"];
    myHoldtime <- BDsummaryStats(mySim)["Holdtime"];
    myStats <- NijBD(mySim);
    if ( identical(myStats, dataCounts) && ##CAREFUL, will this work?identical?
        abs(myHoldtime - dataHoldtime) < eps )
      {return(TRUE);}
    else   {return(FALSE)}
  }
  sapply(sims, ARfn);
}





bdARsimCondEnd.1 <- function(Naccepted=NULL, Ntotal=NULL, Nmax=NULL,
                             T=1.02,
                             L=.3,m=.4,nu=.1, a=8, b=9){
  ##require(acceptRejectSim)
  if (is.null(Naccepted)){
    result <- ARsim(margSimFn=function(){birth.death.simulant(t=T, lambda=L, mu=m, nu=nu, X0=a);},
                    acceptFn=  function(ctmc){list(getMCstate(CTMC=ctmc, Ts=T) == b, NA)},
                    N=Ntotal,keepTestInfo=FALSE);
  }
  else{ #loop until you get at least Naccpted sims. "naccepted*4" arbitrary!
    Ntotal <- 0;
    simsPerRound <- Naccepted*1 ##TOTALLY ARBITRARY; 
    result <- result1 <- vector(mode="list");
    result$acceptSims <- result$testVals <- NULL;
    if (is.null(Nmax) || is.na(Nmax)) Nmax <- Inf
    simCount <- 0;
    while (length(result$acceptSims) < Naccepted &&
           simCount<=Nmax){
      result1 <- ARsim(margSimFn=function(){birth.death.simulant(t=T, lambda=L, mu=m, nu=nu, X0=a);},
                       acceptFn=  function(ctmc){list(getMCstate(CTMC=ctmc, Ts=T) == b, NA)},
                       N=simsPerRound,keepTestInfo=FALSE);
      Ntotal <- Ntotal+simsPerRound; ## may be interested
      result$acceptSims <- c(result$acceptSims, result1$acceptSims);
      simCount <- simCount + simsPerRound
    }
  }
  ##print(Ntotal)
  result$acceptSims; 
}


###THIS IS ACTUALLY NOT TESTED fully.  passes the eyeball test, but didnt
### test this separately from testing bdARsimCondEnd.1 above.

bdARsimCondEnd <- function(Naccepted=NULL, Ntotal=NULL, Nmax=NULL,
                           bd.PO=new("CTMC_PO_1", states=c(5,7,3), times=c(0,.4,1)),
                           L=.5, m=.7, nu=.4){
  if (is(bd.PO, "CTMC_PO_1")){
    numObs <- length(getStates(bd.PO));
    startStates <- getStates(bd.PO)[1:numObs-1]
    endStates <- getStates(bd.PO)[2:numObs]
    startTimes <- getTimes(bd.PO)[1:numObs-1]
    endTimes <- getTimes(bd.PO)[2:numObs]
  }
  else {
    numObs <- length(bd.PO$states);
    startStates <- bd.PO$states[1:numObs-1]
    endStates <- bd.PO$states[2:numObs]
    startTimes <- bd.PO$times[1:numObs-1]
    endTimes <- bd.PO$times[2:numObs]
  }
  deltas <- endTimes-startTimes;
  numInts <- length(deltas);
  simCondArg <- matrix(data=c(startStates, endStates, deltas), ncol=3);
  simCond <- vector("list", numInts);
  simCond <- apply(X=simCondArg, MARGIN= 1,
                   FUN=function(SED){
                     bdARsimCondEnd.1(Naccepted=Naccepted, Ntotal=Ntotal,
                                      Nmax=Nmax,
                                      T=SED[3], a=SED[1], b=SED[2], L=L,
                                      m=m, nu=nu); }
                   );
  if (is.null(Naccepted)){ #so use Ntotal
    #simCond is a list of BDlists; we want the minimum length of BDlist
    numSimConds <- sapply(simCond, length);
    numSimConds <- min(numSimConds);
  }
  else{ #use Naccepted
    numSimConds <- Naccepted
  }
  res <- vector("list", length=numSimConds);
  simCond <- lapply(simCond, function(listLists){listLists[1:numSimConds]}); #truncate
  ##    simCond <- lapply(simCond, function(listLists){listLists$acceptSims[1:numSimConds]}); #truncate
  for (nnn in 1:numSimConds){
    tmp <- vector("list", length=numInts)
    for (iii in 1:numInts){
      tmp[iii] <- (simCond[[iii]])[[nnn]]
    }
    myResult <- combineCTMC(tmp);
    res[nnn] <-
      new("BDMC", times=myResult$times, states=myResult$states, T=myResult$T);
  }
  res;
}




####################################################################
####################################################################
##################################T HIS IS FROM HEMAT CODE
#not sure how it should be organized or where it should go ; 


#pass in BD CTMC; returns times of deaths, times of births
#if getTimes is TRUE returns times as opposed to indices.
#iF its FALSE then return indices (corresponding to $times).
#Returns List with $timesup and $timesdown or $indsup and $indsdown
getBDjTimes <- function(bdMC, getTimes=TRUE){
  #bdMC <- getSubMC(CTMC=bdMC, T=T); #could make this subsample up to T.
  if (is(bdMC, "BDMC")){
    theStates <- getStates(bdMC);
    theTimes <- getTimes(bdMC);
  }
  else{
    theStates <- bdMC$states;
    theTimes <- bdMC$times;
  }
  N <- length(theStates);
  if (N==1){
    if (getTimes){return( list(timesup=numeric(),timesdown=numeric()));}
    else {return( list(indsup=numeric(), indsdown=numeric()) );}
  }#else N>=2
  js <- theStates[2:N] - theStates[1:N-1]; 
  downs <- js == -1;
  inds.down <- seq(1,by=1,to=N-1)[downs] +1;  #+1 to make them indices of $times
  ups <- js == 1;
  inds.up <- seq(1,by=1,to=N-1)[ups] + 1;
  times.down <- theTimes[inds.down];
  times.up <- theTimes[inds.up];
  if(getTimes) {return(list(timesup=times.up,timesdown=times.down));}
  else {return(list(indsup=inds.up, indsdown=inds.down));}
}



#CONDITIONAL ON the first time being up (down)
#ANOTHER way to say it is it either returns the first time or null if that
#first time is in the approrpiate direction
firstJtimeCond <- function(bdMC, up=TRUE){
  theTimes <- getBDjTimes(bdMC, getTimes=TRUE);
  if (up) { #conditional on "goodTimes" being earlier than "badTimes"
    goodTimes <- theTimes$timesup;
    badTimes <- theTimes$timesdown;
  }
  else {
    goodTimes <- theTimes$timesdown;
    badTimes <- theTimes$timesup;
  }
  if ( length(goodTimes)==0 ) { return(NULL);  }
  else if (length(badTimes) == 0) {goodTimes[1]; }
  else if (goodTimes[1] >= badTimes[1]) { return(NULL); }
  else return (goodTimes[1]);
}


