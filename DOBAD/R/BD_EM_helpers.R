
                                        #EM algorithm for birthdeath process.
                                        #"SC"=Special case of nu= n*L
                                        #For the loop to end, either M iterations must be run or lambda and mu
                                        # must be within tol of their previous value.  tol=0 is valid.
                                        #Returns a matrix with the estimator value at each iteration.
                                        #The M+1st index is always the estimator, even if the tolerance stops sooner.
EM.BD.SC.1 <- function(dat, init.params, tol=0.001, M=30, beta.immig, 
                       dr=1e-07, n.fft=1024,r=4,
                       prec.tol=1e-12, prec.fail.stop=TRUE,
                       verbose=1, verbFile=NULL) {
  if (!is.null(verbFile)) sink(verbFile)
  eps <- c(1/0,1/0); #epsilon
  estimators <- init.params;
  estimators.hist <- matrix(nrow=M+1, ncol=2);
  estimators.hist[1,] <- init.params;
  i <- 1;
  while( i <=M && (eps[1]>tol || eps[2]>tol)){
    ##    if (max(estimators*max(dat$times)) > 5) print("Potential problem: estimators are too large.");
    estimators <- getNewParams.SC(oldParams=estimators, beta.immig=beta.immig,theData= dat,
                                  dr=dr, n.fft=n.fft, r=r,prec.tol, prec.fail.stop);
    estimators.hist[i+1,] <- estimators;
    if (verbose>0){
      print( paste("Iteration", i, "just finished and the new estimators are"));
      print(estimators);    
    }
    estimators.hist[M+1,] <- estimators; #for ease of access if tolerance kicks in
    eps <- abs(estimators - estimators.hist[i,]);
    i <- i+1;
  }
  if (!is.null(verbFile)) sink(NULL);
  estimators.hist;
}

getNewParams.SC <- function(theData, oldParams, beta.immig,  dr=0.001,r=4,
                            n.fft=1024, prec.tol, prec.fail.stop){
  UseMethod("getNewParams.SC", theData);
}


E.step.SC.CTMC_PO_many <- function(theData, oldParams, beta.immig,  dr=0.001,
                                   n.fft=1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE){
  n.Procs <- length(theData@BDMCsPO); ##badform to directly reference.
  ##resultMatrix <- sapply(theData@BDMCsPO,
  resultMatrix <- sapply(getBDMCsPOlist(theData), ## this has been modified and not tested
                         function(po1){E.step.SC.CTMC_PO_1(po1,oldParams,beta.immig,
                                                           dr,n.fft,r=r,
                                                           prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)})
  apply(resultMatrix, 1, sum);
}


## This seems to be as fast as the version with 4 calls to mapply, tested for N = 100
## E.step.SC.CTMC_PO_1<- function(theData, oldParams, beta.immig,  dr=0.001,
##                                n.fft=1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE){
##   vec <- c(0,0,0);
##   L <- oldParams[1];
##   mu <- oldParams[2];
##   nu <- beta.immig*L;
##   N <- length(getStates(theData));
##   for (i in 1:(N-1)){
##     timeDiff <- getTimes(theData)[i+1] - getTimes(theData)[i];
##     ##nfftSums <- max(n.fft, getStates(theData)[i+1]+1); #This is quite inefficient.  should code a 1-getter
##     vec[1] <- vec[1]+
##       add.cond.mean.one(lambda=L, mu=mu,
##                         nu=nu, X0=getStates(theData)[i],
##                         t=timeDiff, delta=dr,
##                         Xt=getStates(theData)[i+1], n=n.fft,r=r,
##                         prec.tol=prec.tol, prec.fail.stop=prec.fail.stop);
##     vec[2] <- vec[2] +
##       rem.cond.mean.one(lambda=L, mu=mu,
##                         nu=nu, X0=getStates(theData)[i],
##                         t=timeDiff, delta=dr,
##                         Xt=getStates(theData)[i+1], n=n.fft,
##                         r=r,prec.tol=prec.tol, prec.fail.stop=prec.fail.stop);
##     vec[3] <- vec[3] +
##       timeave.cond.mean.one(lambda=L, mu=mu,
##                             nu=nu, X0=getStates(theData)[i],
##                             t=timeDiff, delta=dr,
##                             Xt=getStates(theData)[i+1], n=n.fft,r=r,
##                             prec.tol=prec.tol, prec.fail.stop=prec.fail.stop);
##   }
##   names(  vec) <- c("Nplus","Nminus", "Holdtime");
##   vec;
## }


## using mapply is not really faster than the for-loop in the old version,
## commented out Have tested that this is identical to old version (with
## for-loop) on a handful of different parameter values, but not
## extraordinarily extensive.  (probably is an all-or-nothign situation
## though)
E.step.SC.CTMC_PO_1<- function(theData, oldParams, beta.immig,  dr=0.001,
                                      n.fft=1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE){
  vec <- c(0,0,0);
  L <- oldParams[1];
  mu <- oldParams[2];
  nu <- beta.immig*L;
  N <- length(getStates(theData));
  arglist <- list(lambda=L,mu=mu,nu=nu,delta=dr,n=n.fft,r=r,prec.tol=prec.tol,
                  prec.fail.stop=prec.fail.stop)
  timeDiff <- diff(getTimes(theData))
  trans.probs <- mapply(process.prob.one, t=timeDiff,
                        X0=getStates(theData)[1:(N-1)], Xt=getStates(theData)[2:N],
                        MoreArgs=list(lambda=L,mu=mu,nu=nu,n=n.fft),
                        SIMPLIFY=TRUE)
  vec[1] <- sum(mapply(add.cond.mean.one, t=timeDiff,X0=getStates(theData)[1:(N-1)], Xt=getStates(theData)[2:N],
                       trans.prob=trans.probs,
                       MoreArgs=arglist, SIMPLIFY=TRUE))
  vec[2] <- sum(mapply(rem.cond.mean.one, t=timeDiff,X0=getStates(theData)[1:(N-1)], Xt=getStates(theData)[2:N],
                       MoreArgs=arglist, SIMPLIFY=TRUE))
  vec[3] <- sum(mapply(timeave.cond.mean.one, t=timeDiff,X0=getStates(theData)[1:(N-1)], Xt=getStates(theData)[2:N],
                       MoreArgs=arglist, SIMPLIFY=TRUE
                       ))
  names(  vec) <- c("Nplus","Nminus", "Holdtime");
  vec;
}


E.step.SC.list <- function(theData, oldParams, beta.immig,  dr=0.001,
                           n.fft=1024,r=4, prec.tol=1e-12, prec.fail.stop=TRUE){
  vec <- c(0,0,0);
  L <- oldParams[1];
  mu <- oldParams[2];
  nu <- beta.immig*L;
  N <- length(getTimes(theData));
  for (i in 1:(N-1)){
    timeDiff <- theData$times[i+1] - theData$times[i];
    nfftSums <- max(n.fft, theData$states[i+1]+1); #This is quite inefficient.  should code a 1-getter
    vec[1] <- vec[1]+
      add.cond.mean.one(lambda=L, mu=mu,
                        nu=nu, X0=theData$states[i],
                        t=timeDiff, delta=dr,
                        Xt=theData$states[i+1], n=n.fft,r=r,
                        prec.tol=prec.tol, prec.fail.stop=prec.fail.stop);
    vec[2] <- vec[2] +
      rem.cond.mean.one(lambda=L, mu=mu,
                        nu=nu, X0=theData$states[i],
                        t=timeDiff, delta=dr,
                        Xt=theData$states[i+1], n=n.fft,r=r
                        ,prec.tol= prec.tol, prec.fail.stop=prec.fail.stop);
    vec[3] <- vec[3] +
      timeave.cond.mean.one(lambda=L, mu=mu,
                            nu=nu, X0=theData$states[i],
                            t=timeDiff, delta=dr,
                            Xt=theData$states[i+1], n=n.fft,r=r,
                            prec.tol=prec.tol, prec.fail.stop=prec.fail.stop);
  }
  names(  vec) <- c("Nplus","Nminus", "Holdtime");
  vec;
}
##Files BD_EM_fns to be not exported in package namespace.
getNewParams.SC.CTMC_PO_many <- function(theData, oldParams, beta.immig,  dr=0.001, n.fft=1024,
                                         r=4,
                                         prec.tol, prec.fail.stop){
  EMsuffStats <- E.step.SC( theData=theData, oldParams=oldParams,
                           beta.immig=beta.immig, dr=dr, n.fft=n.fft,r=r,
                           prec.tol, prec.fail.stop);
  T.total <- sum(sapply(theData@BDMCsPO, function(po1){
    getTimes(po1)[length(getTimes(po1))] - getTimes(po1)[1]}));
  replace(oldParams, 1:2, M.step.SC(EMsuffStats,T.total, beta.immig));
}

getNewParams.SC.CTMC_PO_1<- function(theData, oldParams, beta.immig,  dr=0.001,
                                     r=4,n.fft=1024, prec.tol, prec.fail.stop){
  EMsuffStats <- E.step.SC( theData=theData, oldParams=oldParams, beta.immig=beta.immig, dr=dr, n.fft=n.fft,
                           r=r, prec.tol, prec.fail.stop);
  T <- getTimes(theData)[length(getTimes(theData))] - getTimes(theData)[1]; #latter term should be 0.
  replace(oldParams, 1:2, M.step.SC(EMsuffStats,T, beta.immig));
}

getNewParams.SC.list <- function(theData, oldParams, beta.immig,  dr=0.001, n.fft=1024,
                                 r=4, prec.tol, prec.fail.stop){
  EMsuffStats <- E.step.SC(theData=theData, oldParams=oldParams, beta.immig=beta.immig,  dr=dr, n.fft=n.fft,
                           r=r, prec.tol, prec.fail.stop);
  T <- theData$times[length(theData$times)] - theData$times[1]; #latter term should be 0.
  replace(oldParams, 1:2, M.step.SC(EMsuffStats,T, beta.immig));
}

getNewParams.SC.default <- getNewParams.SC.list

E.step.SC.default <- E.step.SC.list;

## #### This is no faster than with the for loop.  i thought 'for' was slow, but...
## E.step.SC.new2 <- function(oldParams, modelParams, theData, dr=0.001, n.fft=1024){
##   vec <- c(0,0,0);
##   if ( !identical(names(oldParams), c("lambdahat", "muhat")) ){
##     print("Didn't name parameters correctly.");
##     print( names(oldParams) );
##   }
##   L <- oldParams["lambdahat"];
##   mu <- oldParams["muhat"];
##   N <- length(theData$states);
##   timeDiffs <- theData$times[2:N]-theData$times[1:(N-1)];
##   tmp <- apply(matrix(1:(N-1)),1,function(idx){
##     c(add.cond.mean.one(lambda=L, mu=mu,
##                       nu=modelParams["n"]*L, X0=theData$states[idx],
##                       t=timeDiffs[idx], delta=dr,
##                       Xt=theData$states[idx+1], n=n.fft),
##      rem.cond.mean.one(lambda=L, mu=mu,
##                       nu=modelParams["n"]*L, X0=theData$states[idx],
##                       t=timeDiffs[idx], delta=dr,
##                       Xt=theData$states[idx+1], n=n.fft),
##     timeave.cond.mean.one(lambda=L, mu=mu,
##                           nu=modelParams["n"]*L, X0=theData$states[idx],
##                           t=timeDiffs[idx], delta=dr,
##                           Xt=theData$states[idx+1], n=n.fft));});
##   vec <- apply(tmp, 1,sum)
##   names(vec) <- c("Nplus","Nminus", "Holdtime");
##   vec;
## }





getInitParams <-   function(numInitParams=1, summary.PO, T, beta.immig, diffScale){
  initParams <- matrix(data=NA, nrow=6, ncol=2); #6 hardcoded for now ...
  natParams <- M.step.SC(summary.PO, T=T, beta.immig=beta.immig); #"natural"
  names(natParams) <- c("lambdahat", "muhat");
  if (natParams[1] == natParams[2]) {
    natParams[2] <- natParams[2] + diffScale
    natParams[1] <- natParams[1] - diffScale
  } else if (natParams[2] == 0){
    natParams[2] <- .01;
  } else if (natParams[1] == 0){
    natParams[1] = natParams[2] / 2 # do lambda second so its smaller.
  }
  scaleUp <- 3;
  scaleDown <- 1/4
  initParams[1,] <- natParams;
  initParams[4,] <- natParams*scaleUp;
  initParams[3,] <- natParams*scaleDown;
  initParams[2,] <- c(natParams[2],natParams[1]);
  initParams[5,] <- initParams[4,]*scaleUp;
  initParams[6,] <- initParams[4,]*scaleDown;
  initParams <-  initParams[1:numInitParams,];
  if (numInitParams==1) initParams <- t(as.matrix(initParams)); #only way i know
  initParams;
}

