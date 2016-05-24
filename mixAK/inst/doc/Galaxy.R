###
###  Full pdf document describing the code included here is available at
###  http://msekce.karlin.mff.cuni.cz/~komarek/software/mixAK/Galaxy.pdf
###
### ==============================================================================

###################################################
### chunk number 1: What should be created in this Sweave run
###################################################
#line 82 "Galaxy.Rnw"
RUN.TIMECONSUMING.CODE <- FALSE
RUN.ALLOUT <- FALSE


###################################################
### chunk number 2: Directory to store figures
###################################################
#line 90 "Galaxy.Rnw"
FIGDIR <- "./figures/"
FIGKEEPDIR <- "./figuresKeep/"


###################################################
### chunk number 3: Directory with results computed in past
###################################################
#line 101 "Galaxy.Rnw"
RESULTDIR <- "/home/komarek/RESULT_OBJ/mixAK-Galaxy-S081115/"  ### must be changed by the user
RESULT2DIR <- "./RESULT_OBJ/"      ### must be changed by the user


###################################################
### chunk number 4: Options
###################################################
#line 107 "Galaxy.Rnw"
options(width=80)


###################################################
### chunk number 5: Check for existence of directories to store results
###################################################
#line 112 "Galaxy.Rnw"
if (!file.exists(RESULTDIR)){
  stop(paste("Directory ", RESULTDIR, " does not exist.\nYou have to create it or change the value of the variable RESULTDIR.\n"))
}  
if (!file.exists(RESULT2DIR)){
  stop(paste("Directory ", RESULT2DIR, " does not exist.\nYou have to create it or change the value of the variable RESULT2DIR.\n"))
}  


###################################################
### chunk number 6: Load results computed in past
###################################################
#line 122 "Galaxy.Rnw"
if ("Galaxy-Result.RData" %in% dir(RESULT2DIR)){
  load(paste(RESULT2DIR, "Galaxy-Result.RData", sep=""))    
     ## contains RJModel2 (without chains), FixModel2 (without chains), 
     ##          PDensRJ2, PDensFix2  
}else{
  if (!RUN.TIMECONSUMING.CODE){    
    stop(paste("Directory ", RESULT2DIR, " does not contain necessary files.\nSet RUN.TIMECONSUMING.CODE to TRUE.", sep=""))
  }    
}  

if (RUN.ALLOUT){
  if ("Galaxy-RJ2.RData" %in% dir(RESULTDIR)){
    load(paste(RESULTDIR, "Galaxy-RJ2.RData", sep=""))        
       ## contains RJModel2 (contains chains as well)
  }else{  
    if (!RUN.TIMECONSUMING.CODE){    
      stop(paste("Directory ", RESULTDIR, " does not contain necessary files.\nSet RUN.TIMECONSUMING.CODE to TRUE.", sep=""))
    }  
  }          
}


###################################################
### chunk number 7: Load needed packages
###################################################
#line 147 "Galaxy.Rnw"
library("mixAK")
library("coda")


###################################################
### chunk number 8: Read the data and compute a brief summary
###################################################
#line 158 "Galaxy.Rnw"
data("Galaxy", package="mixAK")
summary(Galaxy)


###################################################
### chunk number 9: Draw histogram
###################################################
#line 163 "Galaxy.Rnw"
postscript(paste(FIGDIR, "figGalaxy01.ps", sep=""), width=6, height=6, 
           horizontal=FALSE)
par(mfrow=c(1, 1), bty="n")
hist(Galaxy, prob=TRUE, col="sandybrown", breaks=seq(7, 37, by=0.5), 
     xlab="Velocity (km/sec)", ylab="Density", main="")
dev.off()


###################################################
### chunk number 10: Length of the MCMC
###################################################
#line 184 "Galaxy.Rnw"
nMCMC <- c(burn=100000, keep=500000, thin=10, info=10000)


###################################################
### chunk number 11: Grid of values for the predictive density
###################################################
#line 188 "Galaxy.Rnw"
ygrid <- seq(5, 40, length=500)


###################################################
### chunk number 12: RJ prior distribution
###################################################
#line 201 "Galaxy.Rnw"
RJPrior1 <- list(priorK="uniform", Kmax=30)
if (RUN.TIMECONSUMING.CODE){
  RJModel1 <- NMixMCMC(y0=Galaxy, prior=RJPrior1, nMCMC=nMCMC, 
                       scale=list(shift=0, scale=1), PED=TRUE)
}  


###################################################
### chunk number 13: The prior distribution for the function NMixMCMC was the same as with
###################################################
#line 209 "Galaxy.Rnw"
RJPrior1 <- list(priorK="uniform", Kmax=30,
                 delta=1,
                 priormuQ="independentC", xi=21.7255, D=630.3614,
                 zeta=2, g=0.2, h=0.01586391)


###################################################
### chunk number 14: RJ prior distribution of Richardson and Green
###################################################
#line 222 "Galaxy.Rnw"
RJPrior2 <- list(priorK="uniform", Kmax=30,
                 delta=1,
                 priormuQ="independentC", xi=21.73, D=630.5121,
                 zeta=2*2, g=0.2, h=0.016/2)


###################################################
### chunk number 15: Tuning parameters for RJ-MCMC
###################################################
#line 231 "Galaxy.Rnw"
parRJMCMC2 <- list(par.u1=c(2, 2), 
                   par.u2=c(2, 2), 
                   par.u3=c(1, 1))


###################################################
### chunk number 16: TESTING: RJMCMC
###################################################
#line 241 "Galaxy.Rnw"
set.seed(770328)
TEST.RJModel2 <- NMixMCMC(y0=Galaxy, prior=RJPrior2, RJMCMC=parRJMCMC2, nMCMC=c(burn=100, keep=100, thin=2, info=50), scale=list(shift=0, scale=1), PED=TRUE)
TEST.PDensRJ2 <- NMixPredDensMarg(TEST.RJModel2[[1]], grid=ygrid)  


###################################################
### chunk number 17: RJMCMC
###################################################
#line 249 "Galaxy.Rnw"
if (RUN.TIMECONSUMING.CODE){
  set.seed(770328)
  RJModel2 <- NMixMCMC(y0=Galaxy, prior=RJPrior2, RJMCMC=parRJMCMC2, 
                       nMCMC=nMCMC, scale=list(shift=0, scale=1), PED=TRUE)
}  


###################################################
### chunk number 18: RJMCMC: Acceptance rates of the different move types
###################################################
#line 280 "Galaxy.Rnw"
print(RJModel2[[1]]$moves)
print(RJModel2[[2]]$moves)


###################################################
### chunk number 19: RJMCMC: Basic posterior summary
###################################################
#line 292 "Galaxy.Rnw"
print(RJModel2)


###################################################
### chunk number 20: RJMCMC: Predictive density
###################################################
#line 299 "Galaxy.Rnw"
if (RUN.TIMECONSUMING.CODE){
  PDensRJ2 <- list()
  PDensRJ2[[1]] <- NMixPredDensMarg(RJModel2[[1]], grid=ygrid)
  PDensRJ2[[2]] <- NMixPredDensMarg(RJModel2[[2]], grid=ygrid)
}  


###################################################
### chunk number 21: RJMCMC: Plot of the predictive density
###################################################
#line 309 "Galaxy.Rnw"
postscript(paste(FIGDIR, "figGalaxy02.ps", sep=""), width=6, height=6, horizontal=FALSE)
par(mfrow=c(1, 1), bty="n")
plot(PDensRJ2[[1]], xlab="Velocity (km/sec)")
dev.off()


###################################################
### chunk number 22: RJMCMC: Plot of the conditional predictive densities
###################################################
#line 318 "Galaxy.Rnw"
postscript(paste(FIGDIR, "figGalaxy03.ps", sep=""), width=6, height=6, 
           horizontal=FALSE)
par(mfrow=c(1, 1), bty="n")
plot(PDensRJ2[[1]], K=c(0, 4:7), xlab="Velocity (km/sec)", 
     lty=c(1, rep(2, 4)), col=c("darkblue", rep("red", 4)), 
     lwd=c(2, rep(1, 4)))
dev.off()


###################################################
### chunk number 23: RJMCMC: Plot of the predictive density
###################################################
#line 331 "Galaxy.Rnw"
postscript(paste(FIGDIR, "figGalaxy04.ps", sep=""), width=6, height=10, 
           horizontal=FALSE)
par(mfrow=c(2, 1), bty="n")
#
## Chain 1
hist(Galaxy, prob=TRUE, col="sandybrown", 
     breaks=seq(7, 37, by=0.5), 
     xlab="Velocity (km/sec)", ylab="Density", main="Chain 1")
lines(PDensRJ2[[1]]$x$x1, PDensRJ2[[1]]$dens[[1]], 
      col="darkblue", lwd=2)
#
## Chain 2
hist(Galaxy, prob=TRUE, col="sandybrown", 
     breaks=seq(7, 37, by=0.5), 
     xlab="Velocity (km/sec)", ylab="Density", main="Chain 2")
lines(PDensRJ2[[2]]$x$x1, PDensRJ2[[2]]$dens[[1]], 
      col="darkblue", lwd=2)
dev.off()


###################################################
### chunk number 24: RJMCMC: Choose chain
###################################################
#line 372 "Galaxy.Rnw"
CH <- 1


###################################################
### chunk number 25: RJMCMC: Create mcmc objects
###################################################
#line 377 "Galaxy.Rnw"
if (RUN.ALLOUT){
  start <- RJModel2[[CH]]$nMCMC["burn"] + 1
  end <- RJModel2[[CH]]$nMCMC["burn"] + RJModel2[[CH]]$nMCMC["keep"]
  chK <- mcmc(RJModel2[[CH]]$K, start=start, end=end)
  chgammaInv <- mcmc(RJModel2[[CH]]$gammaInv, start=start, end=end)
  chmixture <- mcmc(RJModel2[[CH]]$mixture, start=start, end=end)
  chdeviance <- mcmc(RJModel2[[CH]]$deviance, start=start, end=end)
}


###################################################
### chunk number 26: RJMCMC: Traceplots
###################################################
#line 392 "Galaxy.Rnw"
if (RUN.ALLOUT){
  lwd <- 0.5
  postscript(paste(FIGKEEPDIR, "figGalaxy07.ps", sep=""), width=6, height=9, 
             horizontal=FALSE)
  par(mfrow=c(2, 2), bty="n")
  traceplot(chK, smooth=FALSE, col="darkgreen", lwd=lwd, main="K")
  traceplot(chgammaInv[, "gammaInv1"], smooth=FALSE, 
            col="brown", lwd=lwd, main="gamma^{-1}")
  traceplot(chmixture[, "y.Mean.1"], smooth=FALSE, 
            col="darkblue", lwd=lwd, main="EY")
  traceplot(chmixture[, "y.SD.1"], smooth=FALSE, 
            col="darkblue", lwd=lwd, main="sd(Y)")
  dev.off()
#
  postscript(paste(FIGKEEPDIR, "figGalaxy08.ps", sep=""), width=6, height=9, 
             horizontal=FALSE)
  par(mfrow=c(2, 2), bty="n")
  traceplot(chdeviance[, "LogL0"], smooth=FALSE, 
            col="red", lwd=lwd, main="Log(L0)")
  traceplot(chdeviance[, "LogL1"], smooth=FALSE, 
            col="red", lwd=lwd, main="Log(L1)")
  traceplot(chdeviance[, "dev.complete"], smooth=FALSE, 
            col="red", lwd=lwd, main="D(complete)")
  traceplot(chdeviance[, "dev.observed"], smooth=FALSE, 
            col="red", lwd=lwd, main="D(observed)")
  dev.off()
}


###################################################
### chunk number 27: RJMCMC: Short traceplot K
###################################################
#line 433 "Galaxy.Rnw"
if (RUN.ALLOUT){
  chKpart <- mcmc(RJModel2[[CH]]$K[490001:500000], start=start+490000, end=end)
  postscript(paste(FIGKEEPDIR, "figGalaxy07a.ps", sep=""), width=9, height=6, 
             horizontal=FALSE)
  par(mfrow=c(1, 1), bty="n")
  traceplot(chKpart, smooth=FALSE, col="darkgreen", main="K")
  dev.off()
}  


###################################################
### chunk number 28: RJMCMC: Density plots
###################################################
#line 455 "Galaxy.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figGalaxy09.ps", sep=""), width=6, height=9, 
             horizontal=FALSE)
  par(mfrow=c(2, 2), bty="n")
  densplot(chK, show.obs=FALSE, col="darkgreen", main="K")
  densplot(chgammaInv[, "gammaInv1"], show.obs=FALSE, 
           col="brown", main="gamma^{-1}", xlim=c(0, 30))
  densplot(chmixture[, "y.Mean.1"], show.obs=FALSE, 
           col="darkblue", main="EY", xlim=c(15, 25))
  densplot(chmixture[, "y.SD.1"], show.obs=FALSE, 
           col="darkblue", main="sd(Y)", xlim=c(0, 12))
  dev.off()
}  


###################################################
### chunk number 29: RJMCMC: Autocorrelation plots
###################################################
#line 478 "Galaxy.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figGalaxy10.ps", sep=""), width=6, height=9, 
             horizontal=FALSE)
  par(mfrow=c(2, 2), bty="n")
  autocorr.plot(chK, auto.layout=FALSE, ask=FALSE, 
                col="darkgreen", lwd=2, main="K")
  autocorr.plot(chgammaInv[, "gammaInv1"], auto.layout=FALSE, ask=FALSE, 
                col="brown", lwd=2, main="gamma^{-1}")
  autocorr.plot(chmixture[, "y.Mean.1"], auto.layout=FALSE, ask=FALSE, 
                col="darkblue", lwd=2, main="EY")
  autocorr.plot(chmixture[, "y.SD.1"], auto.layout=FALSE, ask=FALSE, 
                col="darkblue", lwd=2, main="sd(Y)")
  dev.off()
}  


###################################################
### chunk number 30: Fixed K prior distribution
###################################################
#line 508 "Galaxy.Rnw"
FixPrior2 <- list(priorK="fixed",
                  delta=1,
                  priormuQ="independentC", xi=21.73, D=630.5121,
                  zeta=2*2, g=0.2, h=0.016/2)


###################################################
### chunk number 31: Fixed K MCMC
###################################################
#line 521 "Galaxy.Rnw"
if (RUN.TIMECONSUMING.CODE){
  Keep <- c("iter", "nMCMC", "dim", "prior", "init", "RJMCMC", 
            "scale", "state", "freqK", "propK", "DIC", "moves", 
            "pm.y", "pm.z", "pm.indDev", "pred.dens", "summ.y.Mean", 
            "summ.y.SDCorr", "summ.z.Mean", "summ.z.SDCorr")
  set.seed(770328)  
  FixModel2 <- list()
  PDensFix2 <- list()
  for (k in 1:10){
    cat(paste("K = ", k, "\n-------------------------------\n", sep=""))
    PriorNow <- FixPrior2
    PriorNow$Kmax <- k
    FixModel2[[k]] <- NMixMCMC(y0=Galaxy, prior=PriorNow, nMCMC=nMCMC, 
                               scale=list(shift=0, scale=1), PED=TRUE)
  #
    cat(paste("\nComputation of pred. densities started on ", date(), 
              "\n", sep=""))  
    PDensFix2[[k]] <- list()
    PDensFix2[[k]][[1]] <- NMixPredDensMarg(FixModel2[[k]][[1]], grid=ygrid)
    PDensFix2[[k]][[2]] <- NMixPredDensMarg(FixModel2[[k]][[2]], grid=ygrid)
    cat(paste("Computation of pred. densities finished on ", date(), 
              "\n\n\n", sep=""))
  #
    FixModel2[[k]][[1]] <- FixModel2[[k]][[1]][Keep]
    FixModel2[[k]][[2]] <- FixModel2[[k]][[2]][Keep]  
    class(FixModel2[[k]][[1]]) <- class(FixModel2[[k]][[2]]) <- "NMixMCMC"    
  }
}  


###################################################
### chunk number 32: Fixed K MCMC: Basic posterior summary
###################################################
#line 556 "Galaxy.Rnw"
print(FixModel2[[6]])


###################################################
### chunk number 33: Fixed K: PED and DIC
###################################################
#line 562 "Galaxy.Rnw"
PED <- RJModel2$PED
DIC <- list(Chain1=RJModel2[[1]]$DIC, Chain2=RJModel2[[2]]$DIC)
for (k in 1:length(FixModel2)){
  PED <- rbind(PED, FixModel2[[k]]$PED)
  DIC[[1]] <- rbind(DIC[[1]], FixModel2[[k]][[1]]$DIC)
  DIC[[2]] <- rbind(DIC[[2]], FixModel2[[k]][[2]]$DIC)  
}
rownames(PED) <- rownames(DIC[[1]]) <- rownames(DIC[[2]]) <- c("RJ-MCMC", paste("K = ", 1:length(FixModel2), sep=""))


###################################################
### chunk number 34: Fixed K: print PED
###################################################
#line 573 "Galaxy.Rnw"
print(PED)


###################################################
### chunk number 35: Fixed K: print DIC
###################################################
#line 576 "Galaxy.Rnw"
print(DIC)


###################################################
### chunk number 36: Fixed K MCMC: Plot of the predictive density
###################################################
#line 585 "Galaxy.Rnw"
postscript(paste(FIGDIR, "figGalaxy05.ps", sep=""), width=6, height=6, 
           horizontal=FALSE)
par(mfrow=c(1, 1), bty="n")
hist(Galaxy, prob=TRUE, col="grey90", breaks=seq(7, 37, by=0.5), 
     xlab="Velocity (km/sec)", ylab="Density", main="")
for (k in 1:10){
  lines(PDensFix2[[k]][[1]]$x$x1, PDensFix2[[k]][[1]]$dens[[1]], col="red")
}  
dev.off()


###################################################
### chunk number 37: Fixed K MCMC: Plot of the predictive density
###################################################
#line 599 "Galaxy.Rnw"
postscript(paste(FIGDIR, "figGalaxy06.ps", sep=""), width=6, height=9, 
           horizontal=FALSE)
par(mar=c(3, 2, 2, 1)+0.1)
par(mfrow=c(5, 2), bty="n")
for (k in 1:10){
  hist(Galaxy, prob=TRUE, col="lightblue", breaks=seq(7, 37, by=0.5), 
       xlab="", ylab="", main=paste("K = ", k, sep=""))
  lines(PDensFix2[[k]][[1]]$x$x1, PDensFix2[[k]][[1]]$dens[[1]], col="red", lwd=2)
}  
dev.off()


###################################################
### chunk number 38: Save results
###################################################
#line 632 "Galaxy.Rnw"
if (RUN.TIMECONSUMING.CODE){
  save(list="RJModel2", 
       file=paste(RESULTDIR, "Galaxy-RJ2.RData", sep=""))  
  
  Keep <- c("iter", "nMCMC", "dim", "prior", "init", "RJMCMC", 
            "scale", "state", "freqK", "propK", "DIC", "moves", 
            "pm.y", "pm.z", "pm.indDev", "pred.dens", "summ.y.Mean", 
            "summ.y.SDCorr", "summ.z.Mean", "summ.z.SDCorr")
  #
  RJModel2[[1]] <- RJModel2[[1]][Keep]
  RJModel2[[2]] <- RJModel2[[2]][Keep]
  class(RJModel2[[1]]) <- class(RJModel2[[2]]) <- "NMixMCMC"
  #
  save(list=c("RJModel2", "PDensRJ2", "FixModel2", "PDensFix2"), 
       file=paste(RESULT2DIR, "/Galaxy-Result.RData", sep=""))    
}  


