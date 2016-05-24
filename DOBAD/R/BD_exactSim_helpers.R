# "_helpers" not included in package namespace

p.i <- function(T, a,b,up, L, m, nu, n.fft=1024,
                subdivisions=100){
 integrate(f.i, lower=0, upper=T,
           T=T,a=a,b=b,up=up, L=L,m=m,nu=nu,n.fft=n.fft,
           subdivisions=subdivisions)
##   tryCatch(  integrate(f.i, lower=0, upper=T,
##                        T=T,a=a,b=b,up=up, L=L,m=m,nu=nu,n.fft=n.fft,
##                        subdivisions=subdivisions),
##            error=function(e){
##              sink("DEBUGp.i")
##              print(paste("myU,T,a,b,up,L,m,nu are:",myU,T,a,b,up,L,m,nu))
##              sink(NULL);
##              e;
##              recover();
##            },
##            finally=NULL
##            );
}


#for t as a vector
f.i <- function(t, T, a,b, up, L, m, nu, n.fft=1024){
  if (up==TRUE){i <- a+1; lambda.ai <- a*L + nu}
  else {
    if (a<=0) {return(rep(0, length(t)));} ## not sure if it should error or return 0
    #if (a<=0) {stop(); print("error: a<=0, trying to jumpdown")}
    i <- a-1; lambda.ai <- a*m;
  }
  lambda.a <- a*(L+m) + nu;
  waitProb <- exp(-lambda.a*t); #prob jumping at or after time=t
  pab <- process.prob.one(t=T, lambda=L, mu=m, nu=nu, X0=a, Xt=b, n=n.fft);  
  pib <- process.prob.one(t=T-t, lambda=L, mu=m, nu=nu, X0=i, Xt=b, n=n.fft);
  theDensity <- waitProb * lambda.ai  * pib / pab;
  theDensity[t<0 | t>T] <- rep(0, length=sum(t<0|t>T));
  theDensity;
}





f.i.notNormalized <- function(t, T, a,b, up, L, m, nu, n.fft=1024){
  if (up==TRUE){i <- a+1; lambda.ai <- a*L + nu}
  else {
    if (a<=0) {return(rep(0, length(t)));} ## not sure if it should error or return 0
    #if (a<=0) {stop(); print("error: a<=0, trying to jumpdown")}
    i <- a-1; lambda.ai <- a*m;
  }
  lambda.a <- a*(L+m) + nu;
  waitProb <- exp(-lambda.a*t); #prob jumping at or after time=t
#  pab <- process.prob.one(t=T, lambda=L, mu=m, nu=nu, X0=a, Xt=b, n=n.fft);
  pib <- process.prob.one(t=T-t, lambda=L, mu=m, nu=nu, X0=i, Xt=b, n=n.fft);
  theDensity <- waitProb * lambda.ai  * pib; ### / pab; #"notNormalized"
  theDensity[t<0 | t>T] <- rep(0, length=sum(t<0|t>T));
  theDensity;
  ## "mult" function has yet to make a difference as far as i can tell.
##   theDensity <-
##     apply(matrix(waitProb*lambda.ai*pib),1,
##           function(xx){
##             mult(xx,1, pab, -1); 
##           });
}


######################### helpers for sim.condBD #####

#Simulate a birth-death process conditional on discretely observed data (no "T" wanted)
sim.condBD.main <- function(bd.PO=list(states=c(5,7,3), times=c(0,.4,1)),
                       L=.5, m=.7, nu=.4,  n.fft=1024){
  UseMethod("sim.condBD.main", bd.PO);
}



chooseSim.condBD.1.CTMC_PO_1 <- function(bd.PO=new("CTMC_PO_1", states=c(5,7,3), times=c(0,.4,1)),
                                       L=.5, m=.7, nu=.4,
                                       maxARsims=100,
                                       n.fft=1024){
  arg <- CTMCPO2indepIntervals.CTMC_PO_1(bd.PO)
  simCond <- apply(X=arg, MARGIN=1,
                   FUN=function(SED){chooseSim.condBD.1(T=SED[3],a=SED[1],b=SED[2],L=L,m=m,nu=nu,n.fft=n.fft,maxARsims=maxARsims)})
  res <- combineCTMC(simCond)
  new("BDMC",times=res$times,states=res$states,T=res$T)
}

chooseSim.condBD.CTMC_PO_1 <- function(N=10,
                                       bd.PO=new("CTMC_PO_1", states=c(5,7,3), times=c(0,.4,1)),
                                       L=.5, m=.7, nu=.4,
                                       maxARsims=100,
                                       n.fft=1024){
  return(replicate(n=N, chooseSim.condBD.1.CTMC_PO_1(bd.PO, L, m,nu,maxARsims,n.fft),
                   simplify=FALSE))
}


sim.condBD.main.CTMC_PO_1 <- function(bd.PO=new("CTMC_PO_1", states=c(5,7,3), times=c(0,.4,1)),
                                      L=.5, m=.7, nu=.4,  n.fft=1024){
  ##Replace this with 'CTMCPO 2 indep intervals".
  numObs <- length(getStates(bd.PO));
  startStates <- getStates(bd.PO)[1:(numObs-1)]
  endStates <- getStates(bd.PO)[2:numObs]
  startTimes <- getTimes(bd.PO)[1:(numObs-1)]
  endTimes <- getTimes(bd.PO)[2:numObs]
  deltas <- endTimes-startTimes;
  simCondArg <- matrix(data=c(startStates, endStates, deltas), ncol=3);
  simCond <- apply(X=simCondArg, MARGIN= 1,
                   FUN=function(SED){
                     sim.condBD.1(T=SED[3], a=SED[1], b=SED[2], L=L,
                                  m=m, nu=nu,
                                  n.fft=n.fft)
                   }
                   );
  res <- combineCTMC(simCond);
  new("BDMC", times=res$times, states=res$states, T=res$T)
}

sim.condBD.main.list <- function(bd.PO=list(states=c(5,7,3), times=c(0,.4,1)),
                       L=.5, m=.7, nu=.4,  n.fft=1024){
  numObs <- length(bd.PO$states);
  startStates <- bd.PO$states[1:numObs-1]
  endStates <- bd.PO$states[2:numObs]
  startTimes <- bd.PO$times[1:numObs-1]
  endTimes <- bd.PO$times[2:numObs]
  deltas <- endTimes-startTimes;
  simCondArg <- matrix(data=c(startStates, endStates, deltas), ncol=3);
  simCond <- apply(X=simCondArg, MARGIN= 1,
                   FUN=function(SED){
                     sim.condBD.1(T=SED[3], a=SED[1], b=SED[2], L=L,
                                  m=m, nu=nu,
                                  n.fft=n.fft); }
                   );
  res <- combineCTMC(simCond);
  new("BDMC", times=res$times, states=res$states, T=res$T)  
}

sim.condBD.main.default <- sim.condBD.main.list;

#Simulate a birth-death process conditional on its beginning and ending states
sim.condBD.1 <- function(T=1, a=5, b=7, L=.5, m=.7, nu=.4, n.fft=1024){
  ##print("Exact simulation doesn't work. sample jump time fails.");
  prec <- 1e-9
  BDproc <- list(a,0,T);
  names(  BDproc) <- c("states","times","T");
  timeLeft <- T;
  currState <- a;
  while (timeLeft >0){
    if (currState==b){
      lambda.a <- currState*(L+m) + nu;
      paaT <-  process.prob.one(timeLeft, lambda=L, mu=m, nu=nu, X0=currState, Xt=currState, n=n.fft);
      pa <- exp(-lambda.a*timeLeft) / paaT
      if (pa>1-prec) pa <- 1
      else if (pa<prec) pa <- 0
      ####DEBUG
      if ( is.na(rbinom(1,1,pa))){
        sink("DEBUGsimcondbd1")
        print(paste(currState, L,m,nu, timeLeft, pa, paaT));
        sink(NULL);
      }
      ####end DEBUG
      coinFlip <- rbinom(1, 1, pa)   ##### much better to do all at once than repeatedly
      if ( coinFlip == 1) { #done
        timeLeft=0;  #add nothing to markov chain
      }
      else{ #theres a jump        
        piTup <- p.i(T=timeLeft, a=currState,b=b, up=TRUE, L=L, m=m, nu=nu, n.fft=n.fft)$value
        piTdown <- p.i(T=timeLeft, a=currState,b=b, up=FALSE, L=L, m=m, nu, n.fft=n.fft)$value #dont sum to 1
        if ( runif(n=1,min=0,max=1) < piTup/(piTup+piTdown) ){#jump up. #piTup+piTdown is NOT 1-pa! (need also P(tau>T) )!
          i <- currState+1;
          tau <- sampleJumpTime2(T=timeLeft,a=currState,b=b,up=TRUE,L=L,m=m,nu=nu);
        }
        else { #jump down
          i <- currState-1;
          tau <- sampleJumpTime2(T=timeLeft,a=currState,b=b,up=FALSE,L=L,m=m,nu=nu);
        }
        lastIdx <- length(BDproc$times); #same for states
        BDproc$times <- append(BDproc$times, BDproc$times[lastIdx]+tau);
        BDproc$states <- append(BDproc$states, i);
        timeLeft <- timeLeft-tau
        currState <- i;
      }
    }
    else { # a != b; get jump for certain
      #this code should be merged w/ above, i.e. this is "there's a jump" code like above
      piTup <- p.i(T=timeLeft, a=currState,b=b, up=TRUE, L=L, m=m, nu=nu, n.fft=n.fft)$value
                                        #piTdown <- p.i(T=timeLeft, a=currState,b=b, up=FALSE, L=L, m=m, nu, n.fft=n.fft)$value
      if ( runif(n=1,min=0,max=1) < piTup ){#jump up.
        i <- currState+1;
        tau <- sampleJumpTime2(T=timeLeft,a=currState,b=b,up=TRUE,L=L,m=m,nu=nu);
      }
      else { #jump down
        i <- currState-1;
        tau <- sampleJumpTime2(T=timeLeft,a=currState,b=b,up=FALSE,L=L,m=m,nu=nu);
      }
      lastIdx <- length(BDproc$times); #same for states
      BDproc$times <- append(BDproc$times, BDproc$times[lastIdx]+tau);
      BDproc$states <- append(BDproc$states, i);
      timeLeft <- timeLeft-tau
      currState <- i;
    }
    #currState <- i;
  }#end while loop
  BDproc;
}



sampleJumpTime2 <- function(T, a, b, up=TRUE, L, m, nu){
  myU <- runif(1,min=0,max=1);
  myf.i.unnorm <- function(myt){ f.i.notNormalized(t=myt, T=T, a=a,b=b,up=up,L=L,m=m, nu=nu); }
  normalizationConst <-integrate(myf.i.unnorm, 0, T)$value # "1"
##   normalizationConst <- 
##     tryCatch(integrate(myf.i.unnorm, 0, T)$value, # "1"
##              error=function(e){
##                print(paste("myU,T,a,b,up,L,m,nu are:",myU,T,a,b,up,L,m,nu))
##                e;
##                recover();
##                },
##              finally=NULL
##              );
  myF.i <- function(myttt){
    (integrate(f=myf.i.unnorm, lower=0, upper=myttt)$value / normalizationConst) - myU
  };
  myf.i <- function(myt){myf.i.unnorm(myt)/normalizationConst}
  res <- NRrootfinder(f=myF.i, fprime=myf.i, start= T/2,
                      tol=T /((L+m+nu)*(a+b)* 10000), maxiters=10,min=0,max=T);
  if (res < 0){  ##Note: it's OK for res(ult) to be larger than T.  
    res <- 0
  }
  res;
}
## #Function to sample from the continuous density of the jump time, given
## #whether the jump is up or not.  If a==b then should divide f.i(t) by 1-pa
## ########333TECHNICALLY need to divide by p_i but i think with rejecgtion samplingdont
## #KEEP IN MIND: if you dividew by p_i, add it to all f.i calls!
## sampleJumpTime <- function(T, a, b, up=TRUE, L, m, nu){
##   #one thought: If the max is always at t=0, then this is unnecessary...
##   #get max value of function
## #  N <- 30; #arbitrary
##   multiplier <- 1; #arbitrary. 1.2 is probably way overkill.
##                                         #  aSeq <- seq(from=0,to=T, by=1/N); #arbitrary
##   aSeq=c(0,T);
##   fiVals <- f.i(t=aSeq, T=T, a=a,b=b,up=up,L=L,m=m, nu=nu);
##   max.fi <- max(fiVals); 
##   M <- max.fi*multiplier;
##   tau <- NULL
##   while(is.null(tau)){
##     xSim <- runif(n=1,min=0,max=T);
##     ySim <- runif(n=1,min=0,max=M);
##     if ( ySim < f.i(t=xSim, T=T, a=a,b=b,up=up,
##                     L=L,m=m, nu=nu)  ){# /p_i unneccesary)
##       tau <- xSim;
##     }
##   }
##   tau;
## }

##### dleete this i think?
## # T is in (0, infinity]. Time that we know the exponential is less than.
## ## not being used, although something like this is probably the bets way.
## ## I believe that using min(lambda, mu) works, i havent seen a counterexample
## ##  but have no proof either yet.  
## expInvCond <- function(p, rate, T){
##   condProb <- 1 - exp(-rate * T); #works for T==Inf
##   - log(1-p*condProb) / rate;
## }
