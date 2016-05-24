### R code from vignette source 'BD_MCMC.Rnw'

###################################################
### code chunk number 1: preamble
###################################################
options(continue="+");
##options(warn=2); #warnings become errors


###################################################
### code chunk number 2: setParams
###################################################
library(DOBAD)
######## Generate the "data"
initstate=7;
set.seed(112);
T=5;  L <- .2; mu <- .4;
beta.immig <- .987;
trueParams <- c(L,mu,beta.immig); names(trueParams) <- c("lambda", "mu","beta") #for saving results
dr <- 0.0000000001; #Need |dr| < |L-mu| always o/w get sqrt(negative).
n.fft <- 1024;
delta <- 1;#play with. or make observation intervals distinct 
dat <- birth.death.simulant(t=T, lambda=L, mu=mu, nu=L*beta.immig, X0=initstate);
fullSummary <- BDsummaryStats(dat); fullSummary
#fullSummary <- BDsummaryStats(dat[[1]]); fullSummary
MLEs <- M.step.SC( EMsuffStats=fullSummary, T=T, beta.immig= beta.immig); MLEs
partialData <- getPartialData( seq(0,T,delta), dat);
observedSummary <- BDsummaryStats.PO(partialData); observedSummary;

##Bayesian parameters
L.mean <- 1; M.mean <- 1.1;
aL <- .02;
bL <-  aL / L.mean
aM <- .022;
bM <- aM / M.mean;
print(paste("Variances are", aL/bL^2, "and", aM/bM^2))

N=10
burn=0


###################################################
### code chunk number 3: runMCMC
###################################################

##Rprof(file="mcmc.rprofout")
timer <- system.time(theMCMC <- BD.MCMC.SC(Lguess=L.mean, Mguess=M.mean,
                                           alpha.L=aL, beta.L=bL, # mean 
                                           alpha.M=aM, beta.M=bM, #mean of 
                                           beta.immig=beta.immig,
                                           data= partialData,
                                           burnIn=burn, N=N));
##Rprof(NULL)
#theMCMC
mean(theMCMC[,1]); #lambda
mean(theMCMC[,2]); #mu
L;
mu;

timer;
options(continue=" "); ##undo the setting we changed at top


###################################################
### code chunk number 4: lambdaPlotCode
###################################################
hist(theMCMC[,1], freq=FALSE, breaks=20,
     xlab="Lambda", ylab = "Density",
     main="Posterior of Lambda")
Lmean <- mean(theMCMC[,1])
abline(col="red", v=Lmean)
abline(col="purple", v=L.mean)
#text(col="red", y=-.3, x=Lmean, labels = "L")
x <- seq(from=0,to=1, by=.01);
y <- dgamma(x, shape=aL, rate=bL)
lines(x,y, col="blue")


###################################################
### code chunk number 5: lambdaPlot
###################################################
hist(theMCMC[,1], freq=FALSE, breaks=20,
     xlab="Lambda", ylab = "Density",
     main="Posterior of Lambda")
Lmean <- mean(theMCMC[,1])
abline(col="red", v=Lmean)
abline(col="purple", v=L.mean)
#text(col="red", y=-.3, x=Lmean, labels = "L")
x <- seq(from=0,to=1, by=.01);
y <- dgamma(x, shape=aL, rate=bL)
lines(x,y, col="blue")


###################################################
### code chunk number 6: muPlotCode
###################################################
hist(theMCMC[,2], freq=FALSE, breaks=20,
     xlab="Mu", ylab = "Density",
     main="Posterior of Mu")
Mmean <- mean(theMCMC[,2])
abline(col="red", v=Mmean)
abline(col="purple", v=M.mean)
x <- seq(from=0,to=1, by=.01);
y <- dgamma(x, shape=aM, rate=bM)
lines(x,y, col="blue")


###################################################
### code chunk number 7: muPlot
###################################################
hist(theMCMC[,2], freq=FALSE, breaks=20,
     xlab="Mu", ylab = "Density",
     main="Posterior of Mu")
Mmean <- mean(theMCMC[,2])
abline(col="red", v=Mmean)
abline(col="purple", v=M.mean)
x <- seq(from=0,to=1, by=.01);
y <- dgamma(x, shape=aM, rate=bM)
lines(x,y, col="blue")


