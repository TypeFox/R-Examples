### R code from vignette source 'BD_EM_wImmig.Rnw'

###################################################
### code chunk number 1: preamble
###################################################
options(continue="+");
##options(warn=2); # warnings become errors


###################################################
### code chunk number 2: setParams
###################################################
library(DOBAD)
set.seed(1156);
initstate=4;
T=25;
L <- .3
mu <- .6
beta.immig <- 1.2;
dr <- 0.000001; #Need |dr| < |L-mu| always o/w get sqrt(negative). (numerical differ'n)
n.fft <- 1024;
trueParams <- c(L,mu);
names(trueParams) <- c("lambda", "mu")



###################################################
### code chunk number 3: BD_EM_wImmig.Rnw:204-225
###################################################

##Get the "data"
dat <- birth.death.simulant(t=T, lambda=L, m=mu, nu=L*beta.immig, X0=initstate);
fullSummary <- BDsummaryStats(dat); 
fullSummary
names(fullSummary) <- c("Nplus", "Nminus", "Holdtime");
MLEs.FullyObserved <- M.step.SC( EMsuffStats=fullSummary, T=T, beta.immig= beta.immig); 
#MLEs
###MLE.FullyObserved are NOT the MLE for the EM, but hopefully close as delta --> 0

#delta <- 2
#observetimes <- seq(0,T,delta)
observetimes <- sort(runif(20,min=0,max=T))
partialData <- getPartialData( observetimes, dat);
T <- getTimes(partialData)[length(getTimes(partialData))]

observedSummary <- BDsummaryStats.PO(partialData); observedSummary;

param0 <- c(.8,.9,1.1); names(param0) <- c("lambdahat", "muhat", "nuhat");
#param0



###################################################
### code chunk number 4: runEM
###################################################

## ##############################################################################
## ##########################Do Generic Optimization
## logLike <- function(rates){
##   BDloglikelihood.PO(partialDat=partialData, L=exp(rates[1]), m=exp(rates[2]),
##                      nu=exp(rates[3]), n.fft=1024);
## }
## genericEstimates <- optim(param0, logLike, 
##                           ##method="L-BFGS-B",
##                           ##lower=c(0.0001, 0.0001, .0001), upper=c(100,100,100),
##                           control=list(fnscale=-1))
## print(genericEstimates <- exp(genericEstimates$par))
## print(logLike(log(genericEstimates)))
## ##############################################################################
## ##########################End Generic Optimization

#########RUN EM
iters <- 1
tol <- .0000005;
##myInitParamMat <- rbind(c(0, 1.46,.65),
##                        c(.43, .4, 1.3));
myInitParamMat <- rbind(c(.25,.26,.15));
emOuts <- DOBAD:::EM.BD(dat=partialData, init.params.mat=myInitParamMat, tol=tol, M=iters, 
                        dr=1e-07, n.fft=1024,
                        alpha=.2, beta=.3, fracSimIncr=3, numMCs.i.start=20,
                        outputBestHist=FALSE)

bestparams <- sapply(emOuts, function(emOut){ emOut[[iters+1]]$newParams});
print(bestparams)
#loglikes <- apply( as.matrix(bestparams),2, function(param){logLike(param)});
#print(loglikes);
########### end Run EM


##save.image("BD_EM_wImmigRnw.rsav");


###################################################
### code chunk number 5: conclusion
###################################################
options(continue=" ");


