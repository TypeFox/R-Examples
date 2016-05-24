### R code from vignette source 'BD_EM.Rnw'

###################################################
### code chunk number 1: preamble
###################################################
##options(continue="+");
##options(warn=2); # warnings become errors


###################################################
### code chunk number 2: sourceCode
###################################################
library(DOBAD)


###################################################
### code chunk number 3: setParams
###################################################

## set.seed(1155)
## initstate=4;
## T=8;
##nn <-  10; ## numobservs
## L <- .5
## mu <- .6
## beta.immig <- 1.2;
## dr <- 0.000001; #Need |dr| < |L-mu| always o/w get sqrt(negative). (numerical differ'n)
## n.fft <- 1024;
## trueParams <- c(L,mu);
## names(trueParams) <- c("lambda", "mu")

set.seed(1155)
initstate=4;
T=4;
nn <-  10; ## numobservs
L <- .55
mu <- .6
beta.immig <- 1.2;
dr <- 0.000001; #Need |dr| < |L-mu| always o/w get sqrt(negative). (numerical differ'n)
n.fft <- 1024;
trueParams <- c(L,mu);
names(trueParams) <- c("lambda", "mu")



###################################################
### code chunk number 4: BD_EM.Rnw:214-226
###################################################
##Get the "data"
dat <- birth.death.simulant(t=T, lambda=L, m=mu, nu=L*beta.immig, X0=initstate);
fullSummary <- BDsummaryStats(dat); 
fullSummary
names(fullSummary) <- c("Nplus", "Nminus", "Holdtime");
MLEs.FullyObserved <- M.step.SC( EMsuffStats=fullSummary, T=T, beta.immig= beta.immig); 
#MLEs
###MLE.FullyObserved are NOT the MLE for the EM, but hopefully close as delta-of-observation --> 0

obsTimes <- cumsum(sort(MCMCpack::rdirichlet(1,rep(2,nn))))*T
partialData <- getPartialData(obsTimes, dat);
observedSummary <- BDsummaryStats.PO(partialData); observedSummary;

###################################################
### code chunk number 5: runEM
###################################################
#########RUN EM

iters <- 1;
tol <- .001;

#### You should uncomment the following 2 commands and run -- commented for speed, for now.
## initParamMat <- getInitParams(numInitParams=1, summary.PO=observedSummary, 
##                                       T=T, beta.immig=beta.immig,
##                                       diffScale=100*dr);

## EMtime <- system.time(estimators.hist <-
##                       EM.BD.SC(initParamMat=initParamMat, M=iters, beta.immig=beta.immig,
##                                dat=partialData, dr=dr, n.fft=n.fft, tol=tol)
##                       )[3];
## EMtime;

#### Generic optimization.
## logLike <- function(rates){
##   BDloglikelihood.PO(partialDat=partialData, L=exp(rates[1]), m=exp(rates[2]),
##                      nu=beta.immig*exp(rates[1]), n.fft=1024);
## }
## genericEstimates <- optim(initParamMat, logLike, 
##                           ##method="L-BFGS-B",
##                           ##lower=c(0.0001, 0.0001, .0001), upper=c(100,100,100),
##                           control=list(fnscale=-1))
## print(genericEstimates <- exp(genericEstimates$par))
## print(logLike(log(genericEstimates)))


##Optimal appears to be: c(.4077229, .8642744)
## Run starting from the optimal to get the right values setup for the CIs:
initParamMat <- matrix(c(.41,.86),nrow=1);
names(initParamMat) <- c("lambdahat","muhat")
iters <- 1;
EMtime <- system.time(estimators.hist <-
                      EM.BD.SC(initParamMat=initParamMat, M=iters, beta.immig=beta.immig,
                               dat=partialData, dr=dr, n.fft=n.fft, tol=tol)
                      )[3];
EMtime;

estimators.hist
Lhat <- estimators.hist[iters+1,1]; Lhat
Mhat <- estimators.hist[iters+1,2]; Mhat
MLEs.FullyObserved;
########### end Run EM


###################################################
### code chunk number 6: getCIs
###################################################
IY.a <- getBDinform.PO(partialData, Lhat=Lhat, Mhat=Mhat, 
                         beta.immig=beta.immig, delta=.001)
print(IY.a);
zScr <- 1.96;
Iinv <- solve(IY.a)
Ldist <- sqrt(Iinv[1,1])*zScr
Mdist <- sqrt(Iinv[2,2])*zScr
CI.L <- c(Lhat-Ldist, Lhat+Ldist);CI.L;
CI.M <- c(Mhat-Mdist, Mhat+Mdist);CI.M;



###################################################
### code chunk number 7: conclusion
###################################################
options(continue=" "); ##undo what we set initially


