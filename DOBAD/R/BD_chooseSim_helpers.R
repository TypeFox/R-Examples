# "_helpers" not included in package namespace





### estimate time it takes for a simCond,
exactSimTime1 <- function(a,b, T,  L=.5, m=.7, nu=.4,  dr=.00001, n.fft=1024){
  #compute time for sample jumptime up
  t.each <- system.time(replicate(10,f.i(t=0, T=T, a=a,b=b,up=TRUE, L=L,m=m, nu=nu)))[1]/ 10
  numcalls <- 400 ##estimate NR at ~7 * (1 to 4 groups of 21 calls to f.i)
  ENplus <- add.cond.mean.one(lambda=L, mu=m, nu=nu, X0=a, t=T, delta=dr, Xt=b, n=n.fft);
  ENminus <- rem.cond.mean.one(lambda=L, mu=m, nu=nu, X0=a, t=T, delta=dr, Xt=b, n=n.fft);
  return((ENplus + ENminus ) * numcalls * t.each)
}

ARsimTime1 <- function(a,b, T, L=.5, m=.7, nu=.4,  n.fft=1024){
  p <- process.prob.one(t=T, lambda=L, mu=m, nu=nu, n=n.fft, X0=a, Xt=b);
  ##time scales linearly in T
  ## multiply by 10 to try to get a nonzero number.
  t.each <- system.time(replicate(100,birth.death.simulant(t=T,X0=a,lambda=L, mu=m, nu=nu)))[1]/100
  print(c(p,1/p, t.each))
  return(t.each/p)
}


ARsimTime <- function(bd.PO=list(states=c(5,7,3), times=c(0,.4,1)),
                       L=.5, m=.7, nu=.4,  n.fft=1024){
  ##use ctmcpo2indepintervals here  
  if (is(bd.PO, "CTMC_PO_1")){
    theStates <- getStates(bd.PO);
    theTimes <- getTimes(bd.PO);
  }
  else {
    theStates <- bd.PO$states;
    theTimes <- bd.PO$times;
  }
  numObs <- length(theStates);
  startStates <- theStates[1:numObs-1]
  endStates <- theStates[2:numObs]
  startTimes <- theTimes[1:numObs-1]
  endTimes <- theTimes[2:numObs]
  deltas <- endTimes-startTimes;
  simCondArg <- matrix(data=c(startStates, endStates, deltas), ncol=3);
  pieceTimes <- apply(simCondArg, 1, function(triple){
    ARsimTime1(a=triple[1], b=triple[2],T=triple[3], L=L, m=m, nu=nu, n.fft=n.fft);                                                    
  });
  sum(pieceTimes);  
}


exactSimTime <- function(bd.PO=list(states=c(5,7,3), times=c(0,.4,1)),
                       L=.5, m=.7, nu=.4,  dr=.00001, n.fft=1024){
##use ctmcpo2indepintervals here  
  if (is(bd.PO, "CTMC_PO_1")){
    theStates <- getStates(bd.PO);
    theTimes <- getTimes(bd.PO);
  }
  else {
    theStates <- bd.PO$states;
    theTimes <- bd.PO$times;
  }
  numObs <- length(theStates);
  startStates <- theStates[1:numObs-1]
  endStates <- theStates[2:numObs]
  startTimes <- theTimes[1:numObs-1]
  endTimes <- theTimes[2:numObs]
  deltas <- endTimes-startTimes;
  simCondArg <- matrix(data=c(startStates, endStates, deltas), ncol=3);
  pieceTimes <- apply(simCondArg, 1, function(triple){
    exactSimTime1(a=triple[1], b=triple[2],T=triple[3],
                  L=L, m=m, nu=nu,dr=dr,  n.fft=n.fft);                                                    
  });
  sum(pieceTimes);  
}
