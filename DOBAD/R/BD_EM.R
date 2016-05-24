##### This is my file for functions for the birth death EM.
                                        #Charles Doss


####many of thes functions can be rewritten to use the
#### CTMCPO2indepIntervals fns now


M.step.SC <- function(EMsuffStats, T, beta.immig){
  newParams <- c(-1,-1); #nonsense values
  ##  newParams <- c( mult(EMsuffStats["Nplus"],1, (EMsuffStats["Holdtime"] + beta.immig*T), -1), 
  ##                 mult(EMsuffStats["Nminus"], 1, EMsuffStats["Holdtime"], -1) );
  newParams <- c( EMsuffStats["Nplus"] / (EMsuffStats["Holdtime"] + beta.immig*T) , 
                 EMsuffStats["Nminus"] / EMsuffStats["Holdtime"] );
  names(newParams) <- c("lambdahat", "muhat");
  newParams;
}

E.step.SC <- function(theData, oldParams, beta.immig,  dr=0.001, n.fft=1024,
                      r=4, prec.tol, prec.fail.stop){
  UseMethod("E.step.SC", theData);
}


                                        #log likelihood for a partially observed linear birth death process
                                        #if you have constant interval observation perhaps using fft would be effective...
BDloglikelihood.PO <-  function(partialDat,L,m, nu, 
                                n.fft=1024){
  UseMethod("BDloglikelihood.PO", partialDat);
}

BDloglikelihood.PO.CTMC_PO_many <- function(partialDat,L,m, nu,  n.fft=1024){
  myBDloglike <- function(dat){BDloglikelihood.PO.CTMC_PO_1(dat, L=L,m=m,nu=nu,n.fft=n.fft);}
  return(sum(sapply(partialDat@BDMCsPO, myBDloglike)));     #bad form to directly access component..
}

BDloglikelihood.PO.CTMC_PO_1 <- function(partialDat,L,m, nu,  n.fft=1024){
  numObs <- length(getStates(partialDat));
  startStates <- getStates(partialDat)[1:numObs-1];
  endStates <- getStates(partialDat)[2:numObs];
  timeInts <- getTimes(partialDat)[2:numObs]-getTimes(partialDat)[1:numObs-1]; #will all be the same often
  theArg <- matrix( c(startStates,endStates, timeInts), nrow=numObs-1)
  transProbs <- apply(theArg, 1, function(arg){
    process.prob.one(t = arg[3], lambda = L, mu = m, 
                     nu = nu, X0 = arg[1], Xt = arg[2], n = n.fft);})  
  logLikelihood <- sum(log(transProbs));  
  logLikelihood;
}


BDloglikelihood.PO.list <- function(partialDat,L,m,nu, 
                                    n.fft=1024){
  numObs <- length(partialDat$states);
  startStates <- partialDat$states[1:numObs-1];
  endStates <- partialDat$states[2:numObs];
  timeInts <- partialDat$times[2:numObs]-partialDat$times[1:numObs-1]; #will all be the same often
  theArg <- matrix( c(startStates,endStates, timeInts), nrow=numObs-1)
  transProbs <- apply(theArg, 1, function(arg){
    process.prob.one(t = arg[3], lambda = L, mu = m, 
                     nu = nu, X0 = arg[1], Xt = arg[2], n = n.fft);})
  logLikelihood <- sum(log(transProbs));
  logLikelihood;
}

BDloglikelihood.PO.default <- BDloglikelihood.PO.list;


## see squarem() library(SQUAREM)
EM.BD.SC <- function(dat, initParamMat, tol=1e-4,M=30, beta.immig, 
                     dr=1e-7, n.fft=1024, r=4, 
                     prec.tol=1e-12, prec.fail.stop=TRUE,
                     verbose=1, verbFile=NULL){
  initParamMat2 <- as.data.frame(t(initParamMat))
  estimatorHists <- lapply(initParamMat2,
                           function(initParms){
                             names(initParms) <- c("lambdahat", "muhat"); #necessary ?
                             print(initParms);
                             EM.BD.SC.1(dat=dat, init.params=initParms , beta.immig=beta.immig, M=M, tol=tol,
                                        dr=dr,  r=r,n.fft=n.fft, prec.tol=prec.tol, prec.fail.stop,
                                        verbose=verbose, verbFile=verbFile);
                           });
  estimators <- t(lapply(estimatorHists, function(aMat){aMat[M+1,];}));
  logLikes <- lapply(estimators,  function(params){
    BDloglikelihood.PO(partialDat=dat, L=params[1], m=params[2], nu=params[1]*beta.immig,
                       n.fft=n.fft);});
  maxLike <- which.max(logLikes);
  estimatorHists[[maxLike]];
}







                                        #this is specific case of "Nij" function below rewritten for B-D proc
                                        #just output a 2xmaxval matrix since jumps are only up or down by 1.
                                        # returns a vector, jumpCounts.  jumpCounts[,] is jumps DOWNby1, [2,] is jumps UP.
                                        # BDhist can be av ector of states or a BDMC in class or list form or CTMC_PO_many.
NijBD <- function(BDhist){
  if (is(BDhist, "BDMC")){
    myStates <- getStates(BDhist);
  }
  else if (is(BDhist, "list")){
    myStates <- BDhist$states
  }
  else {#should be vector
    myStates <- BDhist;
  }
  result <- matrix(0, 2, max(myStates)+1);#both dims are same.Excess of 1x2...
  n <- length(myStates);
  currIdx <- myStates[1]+1;
  currS <- myStates[1]
  nextIdx <- NULL;
  nextS <- NULL;  
  for (count in 1:(n-1)){
    nextIdx <- myStates[count+1]+1;
    nextS <- myStates[count+1];
    updownIdx <- (nextS-currS + 3)/2; #map (-1,1) to (1,2)
    result[updownIdx,currIdx] <- result[updownIdx,currIdx]+1;
    currIdx <- nextIdx;
    currS <- nextS;
  }
  result;
}

NijBD.CTMC_many <- function(BDhists){ #should really be BDMC_many not CTMC_many
  ##   if (is(BDhist, "BDMC_many")){
  ##   }
  ##   else if (is(BDhist, "list")){
  ##   }
  maxState <- max(sapply(BDhists@CTMCs, function(ctmc){max(getStates(ctmc))}))
  rowLength <- maxState+1; ## b/c 0 state.
  NijList <- lapply(BDhists@CTMCs, NijBD)
  extendMatrixRows <- function(mat){
    currRowLength <- length(mat[1,]);
    zeroMat <- matrix(data=rep(0,2*(rowLength-currRowLength)),nrow=2,
                      ncol=rowLength-currRowLength)
    cbind(mat,zeroMat)
  }
  NijList <- lapply(NijList, extendMatrixRows)
  NiUps <- sapply(NijList, function(nij){nij[2,]})
  NiDowns <- sapply(NijList, function(nij){nij[1,]})
  NiUp <- apply(NiUps, 1, sum)
  NiDown <- apply(NiDowns,1,sum);
  matrix( c(NiDown,NiUp), nrow=2, byrow=TRUE);
}


                                        #Gets waiting/holding times in each state for a CTMC
                                        #Timehist should have one less entry that stateHist ,ie
                                        # time between stateHist1 and stateHist2 is timeHist1.
                                        #T is the total observed time.
waitTimes <- function(stateHist, timeHist,T){
  result <- vector(mode="numeric", length=max(stateHist)+1);
  n <- length(stateHist);
  currIdx <- stateHist[1]+1;

  count <- 1
  while(count <= n-1){ #works even if n=0
    result[currIdx] <- result[currIdx] + timeHist[count+1] - timeHist[count];
    currIdx <- stateHist[count+1]+1;
    count <- count +1;
  }
  result[currIdx] <- result[currIdx] + T - timeHist[n];
  result;
}


## Note: Nplus, Nminus work with CTMC_PO_1 as well as bdmcs.
## i.e. any class that has a getStates method should work with Nplus.
## Hence, it has no "type" associated.  the "CTMC_PO_many" class
## doesn't have a getStates, so it gets a new function/new name.
Nplus <- function(sim){
  ##sum(NijBD(sim)[2,])
  sum(diff(getStates(sim))>0)
}

Nplus.CTMC_PO_many <- function(ctmcpomany){
  sum(sapply(ctmcpomany@BDMCsPO,Nplus))
}

##sum(sapply(dat@BDMCsPO,Nplus))

## This version seems to work more generally
## Nplus.CTMC_PO_1 <- function(ctmcpo){
##     sum(diff(getStates(ctmcpo))<0)
## }

Nminus <- function(sim){
  ##sum(NijBD(sim)[1,])
  sum(diff(getStates(sim))<0)
}
Nminus.CTMC_PO_many <- function(ctmcpomany){
  sum(sapply(ctmcpomany@BDMCsPO,Nminus))
}


holdTime <- function(sim){
  wts <- waitTimes(getStates(sim),getTimes(sim),getT(sim))
  (0:(length(wts)-1)) %*% wts
}























############################ MCEM
############################


## initparmmat doesnt need names
EM.BD <- function(dat, init.params.mat, tol=0.001, M=30, 
                  dr=1e-07, n.fft=1024,
                  alpha=.2, beta=.3, fracSimIncr=3, numMCs.i.start=50,
                  outputBestHist=TRUE,
                  verbose=1, verbFile=NULL){
  if (!is.null(verbFile)) sink(verbFile)
  initParamMat2 <- as.data.frame(t(init.params.mat))
  EMouts <- lapply(initParamMat2,
                   function(initParms){
                     names(initParms) <- c("lambdahat", "muhat", "nuhat"); #necessary                                                          
                     print(initParms);
                     if ( is(dat,"CTMC_PO_many"))
                       EM.BD.1.CTMC_PO_many(dat, init.params=initParms, tol, M, 
                                            dr, n.fft,
                                            alpha, beta, fracSimIncr, numMCs.i.start,
                                            verbose=verbose)
                     else if ( is(dat,"CTMC_PO_1"))
                       EM.BD.1.CTMC_PO_1(dat, init.params=initParms, tol, M, 
                                         dr, n.fft,
                                         alpha, beta, fracSimIncr, numMCs.i.start,
                                         verbose=verbose)
                   });
  estimators <- t(lapply(EMouts, function(emOut){emOut[[M+1]]$newParams;}));
  logLikes <- lapply(estimators,  function(params){
    BDloglikelihood.PO(partialDat=dat, L=params[1], m=params[2], nu=params[3],
                       n.fft=n.fft);});
  maxLike <- which.max(logLikes);
  if (!is.null(verbFile)) sink(NULL)
  if(outputBestHist) return(EMouts[[maxLike]])
  else return(EMouts)
}


                                        ##initparms needs names!
                                        ##alpha = level/size, beta = type II error = 1-power
EM.BD.1.CTMC_PO_many <- function(dat, init.params, tol=0.001, M=30, 
                                 dr=1e-07, n.fft=1024,
                                 alpha=.2, beta=.3, fracSimIncr=3, numMCs.i.start=50,
                                 verbose=1) {
  ##  eps <- c(1/0,1/0); #epsilon
  names(init.params) <- c("lambdahat", "muhat", "nuhat")
  estimators <- init.params;
  estimators.hist <- matrix(nrow=M+1, ncol=3);
  estimators.hist[1,] <- init.params;
  simOutput.hist <- vector(length=M+1, mode="list");
  simOutput.hist[[1]] <- list(newParams=init.params, deltaQ=NA, sigmaSq=NA);
  NRtol=1e-7; ##will get revised to be related to how much the EM improves
  i <- 1;
  

  zalpha <- qnorm(1-alpha, mean=0,sd=1);
  zbeta <- qnorm(1-beta, mean=0,sd=1);

  upperBound <- 1/0; #probable upper bound for amount of increase in Q function.
  
  ##this is # mcs to star tthe 't'th round with. may need more for the t'th roudn to complete.
                                        #  while( i <=M && (eps[1]>tol || eps[2]>tol)){
  while( i <=M && upperBound > tol ){ #
                                        #    if (max(estimators*max(dat$times)) > 5) print("Potential problem: estimators are too large.");

    if (verbose>0){
      print(paste("The estimators at step", i, "are")); print(estimators);
    }
    simOutput.hist[[i+1]] <- getNewParams.CTMC_PO_many(oldParams=estimators,theData= dat, numMCs.init=numMCs.i.start,
                                          dr=dr, myAlpha=.2, fracSimIncr=fracSimIncr, NRtol=NRtol);
    simOutput.hist[[M+1]] <- simOutput.hist[[i+1]]
    estimators <- simOutput.hist[[i+1]]$newParams;
    tmpSignCheck <- estimators>=0;
    if ( sum(estimators[!tmpSignCheck]) < -1e-8){ ###crude standard; not worth worrying about miniscule float-point errors
      print(paste("EM.BD.1: error in sign of estimators. The estimators are currently", estimators));
    }
    estimators <- tmpSignCheck * estimators; ## replace negatives with 0.        
    deltaQ <- simOutput.hist[[i+1]]$deltaQ
    estimators.hist[i+1,] <- estimators;
###use 'max' for NRtol; note 'min' may yield NRtol==0, which is bad.
    NRtol <- max(abs(estimators.hist[i+1,]- estimators.hist[i,]))/100; #tolerance for newton-raphson maximization;
    estimators.hist[M+1,] <- estimators; #for ease of access if tolerance kicks in
    ##eps <- abs(estimators - estimators.hist[i,]);
    numMCs.i.start <- max(simOutput.hist[[i+1]]$numMCs,
                          simOutput.hist[[i+1]]$sigmaSq*(zalpha+zbeta)^2 / deltaQ^2 )
    upperBound <- deltaQ + zalpha*simOutput.hist[[i+1]]$ASE;
    print(paste("numMCs for next time is ", numMCs.i.start))
    i <- i+1;
  }  
                                        #  estimators.hist;
  simOutput.hist;
}

EM.BD.1.CTMC_PO_1<- function(dat, init.params, tol=0.001, M=30, 
                    dr=1e-07, n.fft=1024,
                    alpha=.2, beta=.3, fracSimIncr=3, numMCs.i.start=50,
                    verbose=1) {
  ##  eps <- c(1/0,1/0); #epsilon
  names(init.params) <- c("lambdahat", "muhat", "nuhat")
  estimators <- init.params;
  estimators.hist <- matrix(nrow=M+1, ncol=3);
  estimators.hist[1,] <- init.params;
  simOutput.hist <- vector(length=M+1, mode="list");
  simOutput.hist[[1]] <- list(newParams=init.params, deltaQ=NA, sigmaSq=NA);
  NRtol=1e-7; ##will get revised to be related to how much the EM improves
  i <- 1;
  

  zalpha <- qnorm(1-alpha, mean=0,sd=1);
  zbeta <- qnorm(1-beta, mean=0,sd=1);

  upperBound <- 1/0; #probable upper bound for amount of increase in Q function.
  
  ##this is # mcs to star tthe 't'th round with. may need more for the t'th roudn to complete.
                                        #  while( i <=M && (eps[1]>tol || eps[2]>tol)){
  while( i <=M && upperBound > tol ){ #
                                        #    if (max(estimators*max(dat$times)) > 5) print("Potential problem: estimators are too large.");
    
    if (verbose>0){
      print(paste("The estimators at step", i, "are")); print(estimators);
    }
    simOutput.hist[[i+1]] <- getNewParams.CTMC_PO_1(oldParams=estimators,theData= dat, numMCs.init=numMCs.i.start,
                                                    dr=dr, myAlpha=.2, fracSimIncr=fracSimIncr, NRtol=NRtol);
    simOutput.hist[[M+1]] <- simOutput.hist[[i+1]]
    estimators <- simOutput.hist[[i+1]]$newParams;
    tmpSignCheck <- estimators>=0;
    if ( sum(estimators[!tmpSignCheck]) < -1e-8){ ###crude standard; not worth worrying about miniscule float-point errors
      print(paste("EM.BD.1: error in sign of estimators. The estimators are currently", estimators));
    }
    estimators <- tmpSignCheck * estimators; ## replace negatives with 0.        
    deltaQ <- simOutput.hist[[i+1]]$deltaQ
    estimators.hist[i+1,] <- estimators;
###use 'max' for NRtol; note 'min' may yield NRtol==0, which is bad.
    NRtol <- max(abs(estimators.hist[i+1,]- estimators.hist[i,]))/100; #tolerance for newton-raphson maximization;
    estimators.hist[M+1,] <- estimators; #for ease of access if tolerance kicks in
    ##eps <- abs(estimators - estimators.hist[i,]);
    numMCs.i.start <- max(simOutput.hist[[i+1]]$numMCs,
                          simOutput.hist[[i+1]]$sigmaSq*(zalpha+zbeta)^2 / deltaQ^2 )
    upperBound <- deltaQ + zalpha*simOutput.hist[[i+1]]$ASE;
    print(paste("numMCs for next time is ", numMCs.i.start))
    i <- i+1;
  }  
                                        #  estimators.hist;
  simOutput.hist;
}






##two vectors of uneven length by extending one with 0s
add.uneven <- function(x,y){
  if (is.vector(x))    return(add.uneven.vector(x,y))
  else if (is.matrix(x)) return(add.uneven.matrix(x,y))
  else if (is.list(x)) return(add.uneven.list(x,y))
}
"%au%" <- add.uneven;

add.uneven.matrix <- function(x,y){
  if (is.null(x)) return(y)
  else if (is.null(y)) return(x)
  ncx <- ncol(x);
  ncy <- ncol(y);
  nrx <- nrow(x);
  nry <- nrow(y);
  if (ncx==ncy){
    maxn <- max(nrx,nry)
    y <- rbind(y, matrix(rep(0,ncx*(maxn-nry)),ncol=ncx));
    x <- rbind(x,matrix(rep(0,ncx*(maxn-nrx)), ncol=ncx));
  }
  else{
    maxn <- max(ncx,ncy)
    y <- cbind(y, matrix(rep(0,nrx*(maxn-ncy)), nrow=nrx))
    x <- cbind(x, matrix(rep(0,nrx*(maxn-ncx)),nrow=nrx));
  }
  x+y;
}
"%aum%" <- add.uneven.matrix

add.uneven.vector <- function(x,y){
  nx <- length(x);
  ny <- length(y);
  maxn <- max(nx,ny);
  x <- c(x,rep(0,maxn-nx))
  y <- c(y,rep(0,maxn-ny))
  return(x+y);
}
"%auv%" <- add.uneven.vector

add.uneven.list <- function(x,y){ ##both lists
  if (is.null(x)) return(y)
  else if (is.null(y)) return(x)
  n <- length(x); #otherwise should be same length.
  if(n!=length(y)) {print("add.uneven.list error list lengths unequal")}
  res <- vector("list", length=n);
  for (i in 1:n){
    res[[i]] <- x[[i]] %au% y[[i]];
  }
  res
}
"%aul%" <- add.uneven.list

getNewParams.CTMC_PO_many <- function(theData, oldParams, numMCs.init,
                         myAlpha=0.2, fracSimIncr=3, dr=0.001, NRtol=1e-7){
  if (!hasArg(debugMode)){
    debugMode <- FALSE
  }

  debugMode <- TRUE
  
  L <- oldParams["lambdahat"];
  mu <- oldParams["muhat"];
  nu <- oldParams["nuhat"];
  newParams <- oldParams;
  T <- getT(theData); #works for either po_1 or po_many
  datLen <- length(theData@BDMCsPO);
  zalpha <- qnorm(1-myAlpha, mean=0,sd=1);
  jumpCountsList <- NULL;
  nisDlist <- NULL;
  nisUlist <- NULL;
  nisU <- NULL
  nisD <- NULL
  myCounts <- NULL;
  myCountsMat <- NULL;
  xisOld <- xisNew <-  NULL;
  ##  numMCs.t <- numMCs.init / (1+1/fracSimIncr);
  numMCs.t <- numMCs.init;
  lowerBound <- -1;
  sims <- NULL; newSims <- NULL;
  Qnew <- Qold <- ASE <- sigmaSq <- deltaQ <- NULL ##confused abou why this is asked for .. some bug i can't get to yet
  numLoops <- 0;
  ##"ascent based"; choose numSims_t
  while ( lowerBound < 0) {
    numLoops <- numLoops+1;
    if (numLoops==1) {numNewMCs <- ceiling(numMCs.init);}
    else {
      numNewMCs <- ceiling(numMCs.t / fracSimIncr);
      numMCs.t <- numMCs.t + numNewMCs;
    }    
    ## E step
    ## this isn't used yet
    ##nisUmaxLengthOld <- length(nisU); nisDmaxLengthOld <- length(nisD);
    myE.step <- function(ctmcpo1){E.step(numNewMCs, theData=ctmcpo1, L,mu,nu)}
    nisUnew <- nisDnew <- myCountsNew <-
      myCountsMatNew <- jumpCountsListNew <- NULL;
    for (i in 1:datLen){
      newres <- myE.step(theData@BDMCsPO[[i]])
      nisUnew <- nisUnew %auv% newres$nisUnew
      nisDnew <- nisDnew %auv% newres$nisDnew
      myCountsMatNew <- myCountsMatNew %aum% newres$myCountsMatNew
      myCountsNew <- myCountsNew %auv% newres$myCountsNew
      jumpCountsListNew <- newres$jumpCountsListNew %aul% jumpCountsListNew
    }
    numOld <- numMCs.t-numNewMCs; #numNewMCs is numnew
    if (numLoops ==1){
      nisU <- nisUnew; 
      nisD <- nisDnew;
      myCounts <- myCountsNew;
    }
    else {
      nisU <- ((numOld*nisU) %auv% (numNewMCs*nisUnew)) / numMCs.t #sum works w/ NULL
      nisD <- ((numOld*nisD) %auv% (numNewMCs*nisDnew)) / numMCs.t;
      myCounts <- (myCounts*numOld + myCountsNew*numNewMCs) / numMCs.t; #sum works w/ NULL
    }
    jumpCountsList <- c(jumpCountsList, jumpCountsListNew); #need later
    myCountsMat <- cbind(myCountsMat, myCountsMatNew);

    
    ##M-step ;
    ## RECODE WITH "M.STEP" FUNCTION!!!!!
    muHat <- myCounts["Nminus"] / myCounts["Holdtime"];
    ##Note we take negative for minimizing
    ## NOTE: using the analytic gradient and hessian apparently slows it down.
    ## This is not clear to me .....
##     nlmArg <- function(nu){
##       res <- -logL.birth.nu.MLE(nu,nisU,myCounts["Nplus"],myCounts["Holdtime"],T=T);
##       attr(res,"gradient") <- -logL.birth.dNu.MLE(nu=nu, nisUp=nisU, holdTime=myCounts["Holdtime"],
##                                                   Nplus=myCounts["Nplus"], T=T);
##       attr(res,"hessian") <- -logL.birth.d2Nu2.MLE(nu=nu, nisUp=nisU, holdTime=myCounts["Holdtime"],
##                                                    Nplus=myCounts["Nplus"], T=T);
##       ###Probably faster if i leave out the grad and hessian ... 
##       res
##     }
##     logld1nu.1arg <-function(aaa){
##       logL.birth.dNu.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
##                          Nplus=myCounts["Nplus"], T=T);}
##     logld2nu.1arg <-function(aaa){
##       logL.birth.d2Nu2.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
##                            Nplus=myCounts["Nplus"], T=T);}
##     ##    nuHat <- NRrootfinder(f=logld1nu.1arg, fprime=logld2nu.1arg, start=oldParams["nuhat"]);
##     nuHat <- solveNuHat (f=logld1nu.1arg, fprime=logld2nu.1arg, start=oldParams["nuhat"], Nplus=myCounts["Nplus"], T=T,
##                         tol=NRtol, maxiters=500); #prefer to use the tolerance.
    
    ##nuHat <- nlm(f=nlmArg,p=oldParams["nuhat"],iterlim=100, gradtol=NRtol)$estimate

##     if (debugMode==TRUE){
##       save(nisU,myCounts,T,oldParams,NRtol,
##            theData, oldParams, numMCs.init,
##            file="getNewParamsNuHatDebug.rsav")
##     }
    
    nuHat <- solveNuHat(nisU, myCounts, T, nuGuess=oldParams["nuhat"],NRtol=NRtol)

    if (nuHat< .04) {
      print("nuhat<.04, saving output; delete this piece of code eventually")
      save(file="BDEMnuHatDebug.rsav", nisD, nisU=nisU, myCounts=myCounts,T=T,oldParams, NRtol,
           theData,L,mu,nu, numNewMCs)
    }
    
    Lhat <- getLhat(nuHat=nuHat, Nplus=myCounts["Nplus"], holdTime=myCounts["Holdtime"],
                    T=T);
    newParams <- replace(newParams, 1:3, c(Lhat,muHat,nuHat));
    xisNew <- logLvec(newParams,jumpCountList=jumpCountsList,myHoldtimeVec=myCountsMat[3,],T=T);
    xisOld <- logLvec(oldParams,jumpCountList=jumpCountsList,myHoldtimeVec=myCountsMat[3,],T=T);
    sigmaSq <- var(xisNew-xisOld); # sigmaSq goes to 0 as lambda(t)-lambda(t-1) does
    ASE <- sqrt(sigmaSq / numMCs.t);
    ##now decide whether or not we improved; simulate until we do
    ##compute deltaQ, ie improvement (or lack thereof)
    jumpCountAvgs <- matrix( c(nisD, nisU), byrow=TRUE, nrow=2);
    Qnew <- Qexact(newParams=newParams, jumpCounts=jumpCountAvgs,
                   holdTime=myCounts["Holdtime"], T=T);
    Qold <- Qexact(newParams=oldParams, jumpCounts=jumpCountAvgs,
                   holdTime=myCounts["Holdtime"], T=T);
    deltaQ <- Qnew-Qold;
    lowerBound <- deltaQ - zalpha*ASE; #we hope it's positive
    if (!is.numeric(lowerBound)) {
      lowerBound <- -1
      print("lowerBound not numeric: automatically continuing. careful of inf loop.");
    }
    if (!is.numeric(Qold) && is.numeric(Qnew)) {
      lowerBound <- 1
      print("Qold not numeric; can't decide if we improved; giving up on ascent property!")
    }
  }
  list(newParams=newParams, deltaQ=deltaQ, sigmaSq=sigmaSq, ASE=ASE, numMCs=numMCs.t)
}


E.step <- function(numNewMCs, theData, L,mu,nu){
  ##newSims <- sim.condBD(N=numNewMCs, bd.PO=theData,
  ##                      L=L, m=mu, nu=nu);
  ##newSims <- bdARsimCondEnd(Naccepted=numNewMCs, Ntotal=NULL, bd.PO=theData,
  ##                          L=L, m=mu, nu=nu);
  newSims <- chooseSim.condBD.CTMC_PO_1(N=numNewMCs,bd.PO=theData,
                                        L=L,m=mu,nu=nu,
                                        maxARsims=100)

  
  jumpCountsListNew <- lapply(newSims, function(mySim){NijBD(mySim)});
  nisDlistNew <- lapply(jumpCountsListNew, function(mat){mat[1,];})
  nisUlistNew <- lapply(jumpCountsListNew, function(mat){mat[2,];})
  myCountsMatNew <- sapply(newSims, function(mySim){BDsummaryStats(mySim)});
  nisUmaxLengthNew <- max(sapply(nisUlistNew,length));
  nisDmaxLengthNew <- max(sapply(nisDlistNew,length));
  nisUnew <- avgNis(nisUlistNew, nisUmaxLengthNew);
  nisDnew <- avgNis(nisDlistNew, nisDmaxLengthNew);
  myCountsNew <- apply(myCountsMatNew,1, mean)
  if (length(newSims) != numNewMCs) {print("EM: E.step: Help,!")}
  return(list(nisUnew=nisUnew, nisDnew=nisDnew, myCountsNew=myCountsNew,
              myCountsMatNew=myCountsMatNew,
              jumpCountsListNew=jumpCountsListNew));
}





##data should be BDMCmany
##For 3 param model with fully observed data
## This gets the MLE via newton-raphson.
## ie. M.step not SC.
M.step <- function(data, nuGuess=1, NRtol=1e-9){
  T <- getT(data)
  data <- data@CTMCs
  jumpCountsList <- lapply(data, function(mySim){NijBD(mySim)});
  ##jumpCountsList <- lapply(data, NijBD); #think this should work...
  ##nisDlistNew <- lapply(jumpCountsList, function(mat){mat[1,];})
  nisUlistNew <- lapply(jumpCountsList, function(mat){mat[2,];})
  nisUmaxLengthNew <- max(sapply(nisUlistNew,length));
  myCountsMatNew <- sapply(data, function(mySim){BDsummaryStats(mySim)});
  nisUnew <- avgNis(nisUlistNew, nisUmaxLengthNew);
  nisUnew <- avgNis(nisUlistNew, nisUmaxLengthNew) * length(data)
  myCounts <- apply(myCountsMatNew,1, sum)  
  nisU <- nisUnew;  
  ##M-step ;
  muHat <- myCounts["Nminus"] / myCounts["Holdtime"];
  logld1nu.1arg <-function(aaa){ logL.birth.dNu.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
                                                    Nplus=myCounts["Nplus"], T=T);}
  logld2nu.1arg <-function(aaa){ logL.birth.d2Nu2.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
                                                      Nplus=myCounts["Nplus"], T=T);}
  nuHat <- solveNuHat(nisU=nisU, myCounts, T=T,nuGuess=nuGuess,NRtol=NRtol)
  Lhat <- getLhat(nuHat=nuHat, Nplus=myCounts["Nplus"], holdTime=myCounts["Holdtime"],
                  T=T);
  estimates <- c(Lhat,muHat,nuHat);
  names(estimates) <- c("lambdahat", "muhat", "nuhat")
  estimates;
}

############## 5/13 changed how this uses solveNuhat
## ##data should be BDMCmany
## ##For 3 param model with fully observed data
## ## This gets the MLE via newton-raphson.
## ## ie. M.step not SC.
## M.step <- function(data, nuGuess=1, NRtol=1e-9){
##   T <- getT(data)
##   data <- data@CTMCs
##   jumpCountsList <- lapply(data, function(mySim){NijBD(mySim)});
##   ##jumpCountsList <- lapply(data, NijBD); #think this should work...
##   ##nisDlistNew <- lapply(jumpCountsList, function(mat){mat[1,];})
##   nisUlistNew <- lapply(jumpCountsList, function(mat){mat[2,];})
##   nisUmaxLengthNew <- max(sapply(nisUlistNew,length));
##   myCountsMatNew <- sapply(data, function(mySim){BDsummaryStats(mySim)});
##   nisUnew <- avgNis(nisUlistNew, nisUmaxLengthNew);
##   nisUnew <- avgNis(nisUlistNew, nisUmaxLengthNew) * length(data)
##   myCounts <- apply(myCountsMatNew,1, sum)  
##   nisU <- nisUnew;  
##   ##M-step ;
##   muHat <- myCounts["Nminus"] / myCounts["Holdtime"];
##   logld1nu.1arg <-function(aaa){ logL.birth.dNu.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
##                                                     Nplus=myCounts["Nplus"], T=T);}
##   logld2nu.1arg <-function(aaa){ logL.birth.d2Nu2.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
##                                                       Nplus=myCounts["Nplus"], T=T);}
##   nuHat <- solveNuHat(f=logld1nu.1arg, fprime=logld2nu.1arg, start=nuGuess, Nplus=myCounts["Nplus"], T=T,
##                       tol=NRtol, maxiters=500); #prefer to use the tolerance.
##   Lhat <- getLhat(nuHat=nuHat, Nplus=myCounts["Nplus"], holdTime=myCounts["Holdtime"],
##                   T=T);
##   estimates <- c(Lhat,muHat,nuHat);
##   names(estimates) <- c("lambdahat", "muhat", "nuhat")
##   estimates;
## }








#######################new version is tested against this one , works.
## #alpha = level/size, beta = type II error = 1-power
## # fracSimIncr (aka "k" in Caffo, Jank, Jones) is geometric rate of increase.
## getNewParams.old <- function(theData, oldParams, numMCs.init,
##                          myAlpha=0.2, fracSimIncr=3, dr=0.001){



##   if (is(theData, "CTMC_PO_1")){
##     theStates <- getStates(theData);
##     theTimes <- getTimes(theData);
##   }
##   else {
##     theStates <- theData$states;
##     theTimes <- theData$times;
##   }
##   T <- theTimes[length(theTimes)] - theTimes[1]; #latter should be 0
##   L <- oldParams["lambdahat"];
##   mu <- oldParams["muhat"];
##   nu <- oldParams["nuhat"];
##   newParams <- oldParams;

##   zalpha <- qnorm(1-myAlpha, mean=0,sd=1);

##   ##MC, E step
##   ##### Should compare whether this is faster or slower than the simcondmethod.
##   ##  ... ie need to write a "simcond" fn that chooses between the two.
##   ##  and in particular, a "decidewhichsimcond" fn so that i can decide (ahead of time)
##   ## here.


## #  numMCs.t <- numMCs.init / (1+1/fracSimIncr);
##   numMCs.t <- numMCs.init;
##   lowerBound <- -1;
##   sims <- NULL;
##   numLoops <- 0; #for debugging/analysis only
##   ##"ascent based"; choose numSims_t
##   while (lowerBound < 0){
##     numLoops <- numLoops+1;
##     print(paste("numloops is", numLoops)); #analyze initial numsim choice
##     if (numLoops==1) {numNewMCs <- ceiling(numMCs.init);}
##     else {
##       numNewMCs <- ceiling(numMCs.t / fracSimIncr);
##       numMCs.t <- numMCs.t + numNewMCs;
##     }

## #####note EXACT SIM DOESNT OWRK NOW
## #### wil rewrite /uncomment once exactsim (ie samplejumptime) works
##                                         #Really not comparable this is more like driver code
## ##     if (exactSimTime(bd.PO=theData, L=L,m=mu, nu=nu, dr=dr) <
## ##         ARsimTime(bd.PO=theData, L=L,m=mu, nu=nu) ){
## ##       newSims <- sim.condBD(N=numNewMCs,bd.PO=theData, L=L, m=mu, nu=nu);
## ##     }
## ##     else {
## ##       newSims <- bdARsimCondEnd(Naccepted=numNewMCs, Ntotal=NULL, bd.PO=theData,
## ##                                 L=L, m=mu, nu=nu);
## ##     }
##     newSims <- bdARsimCondEnd(Naccepted=numNewMCs, Ntotal=NULL, bd.PO=theData,
##                               L=L, m=mu, nu=nu);    
##     sims <- c(sims, newSims);

##     print(length(newSims));
##     print(numNewMCs);

##     ## ENORMOUSLY INEFFICIENT RIGHT NOW ; no need to recompute all these things
##     ## for each round; can easily just compute it for the new ones

##     ## E step
##     jumpCountsList <- lapply(sims, function(mySim){NijBD(mySim)});
##     nisDlist <- lapply(jumpCountsList, function(mat){mat[1,];})
##     nisUlist <- lapply(jumpCountsList, function(mat){mat[2,];})
##     nisUmaxLength <- max(sapply(nisUlist,length));
##     nisDmaxLength <- max(sapply(nisDlist,length));
##     nisU <- avgNis(nisUlist, nisUmaxLength);
##     nisD <- avgNis(nisDlist, nisDmaxLength);
##     myCountsMat <- sapply(sims, function(mySim){BDsummaryStats(mySim)})
##     myCounts <- apply(myCountsMat,1, mean);
##     ##M-step ;
##     muHat <- myCounts["Nminus"] / myCounts["Holdtime"];
##     logld1nu.1arg <-function(aaa){ logL.birth.dNu.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
##                                                       Nplus=myCounts["Nplus"], T=T);}
##     logld2nu.1arg <-function(aaa){ logL.birth.d2Nu2.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
##                                                         Nplus=myCounts["Nplus"], T=T);}
## #    nuHat <- NRrootfinder(f=logld1nu.1arg, fprime=logld2nu.1arg, start=oldParams["nuhat"]);
##     nuHat <- solveNuHat(f=logld1nu.1arg, fprime=logld2nu.1arg, start=oldParams["nuhat"], Nplus=myCounts["Nplus"], T=T);
##     Lhat <- getLhat(nuHat=nuHat, Nplus=myCounts["Nplus"], holdTime=myCounts["Holdtime"],
##                     T=T);
##     newParams <- replace(newParams, 1:3, c(Lhat,muHat,nuHat));
##     xisNew <- logLvec(newParams,jumpCountList=jumpCountsList,myHoldtimeVec=myCountsMat[3,],T=T);
##     xisOld <- logLvec(oldParams,jumpCountList=jumpCountsList,myHoldtimeVec=myCountsMat[3,],T=T);
##     sigmaSq <- var(xisNew-xisOld); # sigmaSq goes to 0 as lambda(t)-lambda(t-1) does 
##     ASE <- sqrt(sigmaSq / numMCs.t);
##     ##now decide whether or not we improved; simulate until we do
##     ##compute deltaQ, ie improvement (or lack thereof)
##     jumpCountAvgs <- matrix( c(nisD, nisU), byrow=TRUE, nrow=2);
##     Qnew <- Qexact(newParams=newParams, jumpCounts=jumpCountAvgs,
##                    holdTime=myCounts["Holdtime"], T=T);
##     Qold <- Qexact(newParams=oldParams, jumpCounts=jumpCountAvgs,
##                    holdTime=myCounts["Holdtime"], T=T);
##     deltaQ <- Qnew-Qold;
##     lowerBound <- deltaQ - zalpha*ASE; #we hope it's positive
##   }
##   list(newParams=newParams, deltaQ=deltaQ, sigmaSq=sigmaSq, ASE=ASE, numMCs=numMCs.t);
## }











##jumpCOunts[1,] should be jumps down, jumpCounts[2,] jumps up.
## jumpCounts and holdTime should be  expecatations (or approximations thereof)
## "exact" sampling, ie not importance sampling
Qexact <- function(newParams, jumpCounts, holdTime,T){
  n <- length(jumpCounts[1,]);
  L <- newParams["lambdahat"]; m <- newParams["muhat"]; nu <- newParams["nuhat"]
  holdTerm <- -(holdTime*(L + m) + T*nu);
  stateSeq <- seq(from=0, to=length(jumpCounts[1,])-1, by=1);
  jumpTerm <- jumpCounts[2,] %*% log(stateSeq*L + nu) +
    jumpCounts[1,2:n] %*% log(stateSeq[2:n]*m); ## don't want log(0)
  jumpTerm + holdTerm;
}

## they are the same function but we think of them differently;
## the "statistics" are expectatoins for Q, and exact values for logL.
logL <- Qexact;



##only point of this function is because a for loop is about as efficient
## as just using 'apply' since you hvae to merge two lists anyway to use apply
##returns vector of loglikelihoods (of complete data)
logLvec <- function(params,jumpCountList, myHoldtimeVec, T){
  n <- length(jumpCountList); #should be same as length(holdtimelist);
  result <- vector(length=n);
  for (i in 1:n){
##    print(paste("loglvec, i " ,i));
##    print(jumpCountList[[i]])
    result[i] <- logL(newParams=params, jumpCounts=jumpCountList[[i]],
                      holdTime=myHoldtimeVec[i], T=T); 
  }
  return(result);
}



getNuHat <- function(Lhat, Nplus, holdTime, T){
  (Nplus - holdTime*Lhat) / T;
}

getLhat <- function(nuHat, Nplus, holdTime, T){
  (Nplus - T*nuHat) / holdTime;
}

## part of likelihood relating only to birth, ie jumps up; mu can be done separately
## nisup is counts of jumps up; nisup[k] is number jumps up in data from state k-1 
logL.birth <- function(newParams,nisUp, holdTime, T){
  L <- newParams["lambdahat"]; nu <- newParams["nuhat"]
  holdTerm <- -(holdTime*L + T*nu);
  stateSeq <- seq(from=0, to=length(nisUp)-1, by=1);
  jumpTerm <- nisUp %*% log(stateSeq*L + nu);
  jumpTerm + holdTerm;
}




                                        # it's ".MLE" because it's evaluated at the MLE, ie at the point where
                                        # we can write nuhat as a function of Lhat.
### this is deriv w.r.t nu, but _not_ deriv along the MLE surface. ie,
### plug in lambda hat as function of nu and the deriv is different.
## logL.birth.dNu <- function(newParams, nisUp, T){
##   L <- newParams["lambdahat"]; nu <- newParams["nuhat"]
##   stateSeq <- seq(from=0, to=length(jumpCounts[1,])-1, by=1);
##   -T + nisUp %*% (1/ (stateSeq*L + nu));
## }

logL.birth.nu.MLE <- function(nu, nisUp, Nplus=sum(nisUp),  holdTime, T){
  L <- getLhat(nuHat=nu, Nplus=Nplus, holdTime=holdTime, T=T);
  holdTerm <- -(holdTime*L + T*nu);
  stateSeq <- seq(from=0, to=length(nisUp)-1, by=1);
  jumpTerm <- nisUp %*% log(stateSeq*L + nu);
  jumpTerm + holdTerm;
}


## derivative of (likelihood with lambda as a function of nu)
## (ie not just deriv likelihood, and _then_ plug in nu.)
logL.birth.dNu.MLE <- function(nu, nisUp, Nplus=sum(nisUp), holdTime, T){
  L <-   getLhat(nuHat=nu, Nplus=Nplus, holdTime=holdTime,T=T);
  stateSeq <- seq(from=0, to=length(nisUp)-1, by=1);
  nisUp %*% (  (1- stateSeq*T/holdTime)/(stateSeq*L + nu) );  
}

## second deriv "d2Nu2" = second deriv wrt nu
logL.birth.d2Nu2.MLE <- function(nu, nisUp, Nplus=sum(nisUp), holdTime, T){
  L <-   getLhat(nuHat=nu, Nplus=Nplus, holdTime=holdTime,T=T);
  stateSeq <- seq(from=0, to=length(nisUp)-1, by=1);
  - nisUp %*% (  (1- stateSeq*T/holdTime) /(stateSeq*L + nu) )^2;  
}

## E.step <- function(theData, oldParams, dr=0.001, numMCs){
##   vec <- c(0,0,0);
##   if ( !identical(names(oldParams), c("lambdahat", "muhat", "nuhat")) ){
##     print("Didn't name parameters correctly.");
##     print( names(oldParams) );
##   }

##   if (is(theData, "CTMC_PO_1")){
##     theStates <- getStates(theData);
##     theTimes <- getTimes(theData);
##   }
##   else {
##     theStates <- theData$states;
##     theTimes <- theData$times;
##   }
##   L <- oldParams["lambdahat"];
##   mu <- oldParams["muhat"];
##   nu <- oldParams["nuhat"];
##   N <- length(theStates);
  
  

## ####irrelevant keeping around for rubric
##   ##   for (i in 1:(N-1)){
##   ##     timeDiff <- theTimes[i+1] - theTimes[i];
##   ##     nfftSums <- max(n.fft, theStates[i+1]+1); #This is quite inefficient.  should code a 1-getter
##   ##     vec[1] <- vec[1]+
##   ##        add.cond.mean.one(lambda=L, mu=mu,
##   ##                        nu=nu, X0=theStates[i],
##   ##                        t=timeDiff, delta=dr,
##   ##                        Xt=theStates[i+1], n=n.fft);
##   ##     vec[2] <- vec[2] +
##   ##      rem.cond.mean.one(lambda=L, mu=mu,
##   ##                       nu=nu, X0=theStates[i],
##   ##                       t=timeDiff, delta=dr,
##   ##                       Xt=theStates[i+1], n=n.fft);
##   ##    vec[3] <- vec[3] +
##   ##      timeave.cond.mean.one(lambda=L, mu=mu,
##   ##                       nu=nu, X0=theStates[i],
##   ##                       t=timeDiff, delta=dr,
##   ##                       Xt=theStates[i+1], n=n.fft);
##   ##   }
##   ##   names(  vec) <- c("Nplus","Nminus", "Holdtime");
  
##   vec;
## }


##Calls NRrootfinder, ensuring nu stays in the correct domain
## solveNuHat <- function(f, fprime, start, tol=1e-7, maxiters=50,
##                        Nplus, T){
##                                         # note that the NR formula will evaluate for nu < 0;  
##   nu.root <- NRrootfinder(f,fprime,start,tol,maxiters);
##                                         #for lambdahat to be positive, nu must be in [0,N/T]:
##                                         #If there's a root, then concavity implies we move to the nearest endpoint of the domain.
##   if (nu.root<0) return(0)
##   else if (nu.root>Nplus/T) return(Nplus/T)
##   else return(nu.root)
## }

solveNuHat <- function(nisU, myCounts,T, nuGuess, NRtol){
  nlmArg <- function(lnu){
    res <- -logL.birth.nu.MLE(exp(lnu),
                              nisU,myCounts["Nplus"],myCounts["Holdtime"],T=T);    
    attr(res,"gradient") <- -logL.birth.dNu.MLE(nu=exp(lnu),
                                                nisUp=nisU,
                                                holdTime=myCounts["Holdtime"],
                                                Nplus=myCounts["Nplus"], T=T) * exp(lnu);
    attr(res,"hessian") <- -logL.birth.d2Nu2.MLE(nu=exp(lnu),
                                                 nisUp=nisU, holdTime=myCounts["Holdtime"],
                                                 Nplus=myCounts["Nplus"], T=T) * exp(2*lnu) +
                                                   -logL.birth.dNu.MLE(nu=exp(lnu),
                                                                       nisUp=nisU,
                                                                       holdTime=myCounts["Holdtime"],
                                                                       Nplus=myCounts["Nplus"], T=T) * exp(lnu);
### using grad & hessian is slower than hessian but not by much
### as long as there are no overflows.  hence stepmax.  value chosen arbitrarily.
    res
  }
  lognuHat <- nlm(f=nlmArg,p=log(nuGuess) ,iterlim=100, gradtol=NRtol,
                  stepmax=.5)$estimate #### HERE HERE HERE ".5" needs to be modified 
  nuHat <- exp(lognuHat)
  if (nuHat > myCounts["Nplus"]/T)  nuHat <- myCounts["Nplus"]/T
  ##else if (nuHat<0) nuHat <- 0
  return(nuHat)
}

                                        #newton raphson
                                        # using "fact" that it will have zero at nu > 0 .
                                        # our f will be monotonic but no guarantee of convergence.
NRrootfinder <- function(f, fprime, start, tol=1e-7, maxiters=15, min=-Inf,max=Inf){
  ##hist <- vector("numeric",length=maxiters+1);
  ##hist[1] <- start;
  xprev <- start;
  for (i in 1:maxiters){
    xnext <- xprev - f(xprev) / fprime(xprev);
    ##  hist[i+1] <- xnext;
    if (abs(xnext-xprev) < tol){
      ##return(xnext)
      break;
    }
    else if (xnext<min) {xnext <- min} ##infinite loop possibly.
    else if (xnext>max) {xnext <- max} ##infinite loop possibly.
    xprev <- xnext;
  }
  ##if (xnext<0) return(0)
  ##else return(xnext)
  ##print(paste("NRrootfinder:iters is ",i))
  ##  print(hist);
  xnext;
}


avgNis <- function(nisList, maxlength){
  nisExtendedList <- sapply(nisList, function(nis){
    c(nis, rep(0, maxlength-length(nis)))});
  nisExtendedList <- matrix(nisExtendedList, nrow=maxlength)
  ##note: "as.matrix" fails; issue is if nrows is 1 we don't want a vector.
  apply(nisExtendedList, 1, mean)
}






########################### SHOULD RECODE W OBSERVED, ie use "M.step (not SC)
                                        #alpha = level/size, beta = type II error = 1-power
                                        # fracSimIncr (aka "k" in Caffo, Jank, Jones) is geometric rate of increase.
getNewParams.CTMC_PO_1 <- function(theData, oldParams, numMCs.init,
                         myAlpha=0.2, fracSimIncr=3, dr=0.001, NRtol=1e-7){

  if (is(theData, "CTMC_PO_1")){
    theStates <- getStates(theData);
    theTimes <- getTimes(theData);
  }
  else {
    theStates <- theData$states;
    theTimes <- theData$times;
  }
  T <- theTimes[length(theTimes)] - theTimes[1]; #latter should be 0
  L <- oldParams["lambdahat"];
  mu <- oldParams["muhat"];
  nu <- oldParams["nuhat"];
  newParams <- oldParams;

  zalpha <- qnorm(1-myAlpha, mean=0,sd=1);

  ##MC, E step
##### Should compare whether this is faster or slower than the simcondmethod.
  ##  ... ie need to write a "simcond" fn that chooses between the two.
  ##  and in particular, a "decidewhichsimcond" fn so that i can decide (ahead of time)
  ## here.

  jumpCountsList <- NULL;
  nisDlist <- NULL;
  nisUlist <- NULL;
  nisU <- NULL
  nisD <- NULL
  myCounts <- NULL;
  myCountsMat <- NULL;
  

                                        #  numMCs.t <- numMCs.init / (1+1/fracSimIncr);
  numMCs.t <- numMCs.init;
  lowerBound <- -1;
  sims <- NULL; newSims <- NULL;
  numLoops <- 0;
  ##"ascent based"; choose numSims_t
  while (lowerBound < 0){
    numLoops <- numLoops+1;
    if (numLoops==1) {numNewMCs <- ceiling(numMCs.init);}
    else {
      numNewMCs <- ceiling(numMCs.t / fracSimIncr);
      numMCs.t <- numMCs.t + numNewMCs;
    }

    
    ## E step
    newSims <- chooseSim.condBD.CTMC_PO_1(N=numNewMCs, bd.PO=theData,
                                          L,  mu ,nu,maxARsims=100)##100 up for debate!
    ##     newSims <- bdARsimCondEnd(Naccepted=numNewMCs, Ntotal=NULL, bd.PO=theData,
    ##                               L=L, m=mu, nu=nu);
    ##newSims <- sim.condBD(N=numNewMCs, bd.PO=theData,
    ##                      L=L, m=mu, nu=nu);


    ##calc for new sims
    jumpCountsListNew <- lapply(newSims, function(mySim){NijBD(mySim)});
    nisDlistNew <- lapply(jumpCountsListNew, function(mat){mat[1,];})
    nisUlistNew <- lapply(jumpCountsListNew, function(mat){mat[2,];})
    myCountsMatNew <- sapply(newSims, function(mySim){BDsummaryStats(mySim)});
    nisUmaxLengthNew <- max(sapply(nisUlistNew,length));
    nisDmaxLengthNew <- max(sapply(nisDlistNew,length));
    nisUnew <- avgNis(nisUlistNew, nisUmaxLengthNew);
    nisDnew <- avgNis(nisDlistNew, nisDmaxLengthNew);
    myCountsNew <- apply(myCountsMatNew,1, mean)
    numOld <- numMCs.t-numNewMCs; #numNewMCs is numnew
    if (length(newSims) != numNewMCs) {print("Help,!")}
    nisUmaxLengthOld <- length(nisU);
    nisDmaxLengthOld <- length(nisD);
    newUlength <- max(nisUmaxLengthOld, nisUmaxLengthNew);
    newDlength <- max(nisDmaxLengthOld, nisDmaxLengthNew);
    nisU <- c(nisU, rep(0, newUlength-nisUmaxLengthOld));
    nisD <- c(nisD, rep(0, newDlength-nisDmaxLengthOld));
    nisUnew <- c(nisUnew, rep(0, newUlength-nisUmaxLengthNew));
    nisDnew <- c(nisDnew, rep(0, newDlength-nisDmaxLengthNew));


    
    if (numLoops ==1) {
      nisU <- nisUnew;
      nisD <- nisDnew;
      myCounts <- myCountsNew;
    }
    else {
      nisU <- (numOld*nisU+ numNewMCs*nisUnew) / numMCs.t #sum works w/ NULL
      nisD <- (numOld*nisD+ numNewMCs*nisDnew) / numMCs.t;
      myCounts <- (myCounts*numOld + myCountsNew*numNewMCs) / numMCs.t; #sum works w/ NULL
    }
    jumpCountsList <- c(jumpCountsList, jumpCountsListNew); #need later
    myCountsMat <- cbind(myCountsMat, myCountsMatNew);

    
    ##M-step ;
    ## RECODE WITH "M.STEP" FUNCTION!!!!!
    muHat <- myCounts["Nminus"] / myCounts["Holdtime"];
    logld1nu.1arg <-function(aaa){ logL.birth.dNu.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
                                                      Nplus=myCounts["Nplus"], T=T);}
    logld2nu.1arg <-function(aaa){ logL.birth.d2Nu2.MLE(nu=aaa, nisUp=nisU, holdTime=myCounts["Holdtime"],
                                                        Nplus=myCounts["Nplus"], T=T);}
                                        #    nuHat <- NRrootfinder(f=logld1nu.1arg, fprime=logld2nu.1arg, start=oldParams["nuhat"]);

    ########## HERE Changed things without testing; commented this 'solveNuHat" call
    #### and replaced with the one below.  Not sure how this wasn't broken though since
    #### i dont see a definition of solvenuhat that works for the following call.
##     nuHat <- solveNuHat(f=logld1nu.1arg, fprime=logld2nu.1arg, start=oldParams["nuhat"], Nplus=myCounts["Nplus"], T=T,
##                         tol=NRtol, maxiters=500); #prefer to use the tolerance.


    nuHat <- solveNuHat(nisU=nisU, myCounts, T=T,nuGuess=oldParams["nuhat"],NRtol=NRtol)
    
    Lhat <- getLhat(nuHat=nuHat, Nplus=myCounts["Nplus"], holdTime=myCounts["Holdtime"],
                    T=T);
    newParams <- replace(newParams, 1:3, c(Lhat,muHat,nuHat));
    xisNew <- logLvec(newParams,jumpCountList=jumpCountsList,myHoldtimeVec=myCountsMat[3,],T=T);
    xisOld <- logLvec(oldParams,jumpCountList=jumpCountsList,myHoldtimeVec=myCountsMat[3,],T=T);
    sigmaSq <- var(xisNew-xisOld); # sigmaSq goes to 0 as lambda(t)-lambda(t-1) does 
    ASE <- sqrt(sigmaSq / numMCs.t);
    ##now decide whether or not we improved; simulate until we do
    ##compute deltaQ, ie improvement (or lack thereof)
    jumpCountAvgs <- matrix( c(nisD, nisU), byrow=TRUE, nrow=2);
    Qnew <- Qexact(newParams=newParams, jumpCounts=jumpCountAvgs,
                   holdTime=myCounts["Holdtime"], T=T);
    Qold <- Qexact(newParams=oldParams, jumpCounts=jumpCountAvgs,
                   holdTime=myCounts["Holdtime"], T=T);
    deltaQ <- Qnew-Qold;
    lowerBound <- deltaQ - zalpha*ASE; #we hope it's positive
  }
  list(newParams=newParams, deltaQ=deltaQ, sigmaSq=sigmaSq, ASE=ASE, numMCs=numMCs.t);
}




