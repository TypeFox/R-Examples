###
###  Full pdf document describing the code included here is available at
###  http://msekce.karlin.mff.cuni.cz/~komarek/software/mixAK/Faithful.pdf
###
### ==============================================================================

### R code from vignette source 'Faithful.Rnw'

###################################################
### code chunk number 1: What should be created in this Sweave run
###################################################
RUN.TIMECONSUMING.CODE <- FALSE
RUN.ALLOUT <- FALSE


###################################################
### code chunk number 2: Directory to store figures
###################################################
FIGDIR <- "./figures/"
FIGKEEPDIR <- "./figuresKeep/"


###################################################
### code chunk number 3: Directory with results computed in past
###################################################
RESULTDIR <- "/home/komarek/RESULT_OBJ/mixAK-Faithful-S081115/"   ### must be changed by the user
RESULT2DIR <- "./RESULT_OBJ/"      ### must be changed by the user


###################################################
### code chunk number 4: Options
###################################################
options(width=80)


###################################################
### code chunk number 5: Check for existence of directories to store results
###################################################
if (!file.exists(RESULTDIR)){
  stop(paste("Directory ", RESULTDIR, " does not exist.\nYou have to create it or change the value of the variable RESULTDIR.\n"))
}  
if (!file.exists(RESULT2DIR)){
  stop(paste("Directory ", RESULT2DIR, " does not exist.\nYou have to create it or change the value of the variable RESULT2DIR.\n"))
}  


###################################################
### code chunk number 6: Load results computed in past
###################################################
Kshow <- 3

if ("Faithful-Result.RData" %in% dir(RESULT2DIR)){
  load(paste(RESULT2DIR, "Faithful-Result.RData", sep=""))  
     ## contains ModelK (without chains), MPDensModelK, JPDensModelK
  Model0 <- ModelK[[Kshow]]
  MPDensModel0 <- MPDensModelK[[Kshow]]
  JPDensModel0 <- JPDensModelK[[Kshow]]
}else{
  if (!RUN.TIMECONSUMING.CODE){    
    stop(paste("Directory ", RESULT2DIR, " does not contain necessary files.\nSet RUN.TIMECONSUMING.CODE to TRUE.", sep=""))
  }      
}  

if (RUN.ALLOUT){  
  if (paste("Faithful-Model0", Kshow, ".RData", sep="") %in% dir(RESULTDIR)){
    load(paste(RESULTDIR, "Faithful-Model0", Kshow, ".RData", sep=""))  
       ## contains Model0=ModelK[[Kshow]] (chains included)
  }else{
    if (!RUN.TIMECONSUMING.CODE){    
      stop(paste("Directory ", RESULTDIR, " does not contain necessary files.\nSet RUN.TIMECONSUMING.CODE to TRUE.", sep=""))
    }      
  }  
}


###################################################
### code chunk number 7: Load needed packages
###################################################
library("mixAK")
library("coda")
library("colorspace")


###################################################
### code chunk number 8: Read the data and compute a brief summary
###################################################
data("Faithful", package="mixAK")
summary(Faithful)


###################################################
### code chunk number 9: Draw scatterplot and histogram
###################################################
postscript(paste(FIGDIR, "figFaithful01.ps", sep=""), width=7, height=10, 
           horizontal=FALSE)
par(bty="n")
layout(matrix(c(0,1,1,1,1,0, 2,2,2,3,3,3), nrow=2, byrow=TRUE))
plot(Faithful, col="red", pch=16, 
     xlab="Eruptions (min)", ylab="Waiting (min)")
hist(Faithful$eruptions, prob=TRUE, col="sandybrown", 
     xlab="Eruptions (min)", ylab="Density", main="", 
     breaks=seq(1.4, 5.6, by=0.3))
hist(Faithful$waiting, prob=TRUE, col="sandybrown", 
     xlab="Waiting (min)", ylab="Density", main="")
dev.off()


###################################################
### code chunk number 10: Length of the MCMC
###################################################
nMCMC <- c(burn=100000, keep=500000, thin=10, info=10000)


###################################################
### code chunk number 11: Grid of values for the predictive density
###################################################
ygrid <- list(eruptions=seq(1, 6, length=100), 
              waiting=seq(40, 100, length=100))


###################################################
### code chunk number 12: TESTING: Prior distribution and model with three components and posterior predictive density
###################################################
set.seed(770328)
Prior0 <- list(priorK="fixed", Kmax=3)
TEST.Model0 <- NMixMCMC(y0=Faithful, prior=Prior0, nMCMC=c(burn=100, keep=100, thin=2, info=50))
TEST.MPDensModel0 <- NMixPredDensMarg(TEST.Model0[[1]], grid=ygrid)
TEST.JPDensModel0 <- NMixPredDensJoint2(TEST.Model0[[1]], grid=ygrid)


###################################################
### code chunk number 13: Prior distribution
###################################################
Prior0 <- list(priorK="fixed", Kmax=3)


###################################################
### code chunk number 14: Model with three components
###################################################
if (RUN.TIMECONSUMING.CODE){
  set.seed(777988621)  
  Model0 <- NMixMCMC(y0=Faithful, prior=Prior0, nMCMC=nMCMC, PED=TRUE)
}  


###################################################
### code chunk number 15: The prior distribution was the same as with (eval = FALSE)
###################################################
## Prior0 <- list(priorK="fixed", Kmax=3, 
##                delta=1, 
##                priormuQ="independentC",
##                xi=c(-0.1207, -0.1028), D=diag(c(9.4033, 15.1983)),
##                zeta=3, g=0.2, h=c(1.0635, 0.6580))


###################################################
### code chunk number 16: Model0: Basic posterior summary
###################################################
print(Model0)


###################################################
### code chunk number 17: Model0: Marginal predictive densities
###################################################
if (RUN.TIMECONSUMING.CODE){
  MPDensModel0 <- list()
  MPDensModel0[[1]] <- NMixPredDensMarg(Model0[[1]], grid=ygrid)  
  MPDensModel0[[2]] <- NMixPredDensMarg(Model0[[2]], grid=ygrid)  
}  


###################################################
### code chunk number 18: Model0: Plot of the marginal predictive density
###################################################
postscript(paste(FIGDIR, "figFaithful02.ps", sep=""), width=7, height=5, 
           horizontal=FALSE)
plot(MPDensModel0[[1]])
dev.off()


###################################################
### code chunk number 19: Model0: Joint predictive density
###################################################
if (RUN.TIMECONSUMING.CODE){
  JPDensModel0 <- list()
  JPDensModel0[[1]] <- NMixPredDensJoint2(Model0[[1]], grid=ygrid)  
  JPDensModel0[[2]] <- NMixPredDensJoint2(Model0[[2]], grid=ygrid)  
}  


###################################################
### code chunk number 20: Model0: Plot of the joint predictive density
###################################################
postscript(paste(FIGDIR, "figFaithful03.ps", sep=""), width=7, height=5, 
           horizontal=FALSE)
plot(JPDensModel0[[1]])
dev.off()


###################################################
### code chunk number 21: Model0: Plot of the marginal predictive density with the histogram
###################################################
postscript(paste(FIGDIR, "figFaithful04.ps", sep=""), width=7, height=10, 
           horizontal=FALSE)
par(mfrow=c(2, 2), bty="n")
ylimE <- c(0, max(c(MPDensModel0[[1]]$dens[[1]], MPDensModel0[[2]]$dens[[1]])))
ylimW <- c(0, max(c(MPDensModel0[[1]]$dens[[2]], MPDensModel0[[2]]$dens[[2]])))
for (CH in 1:2){
  hist(Faithful$eruptions, prob=TRUE, col="sandybrown", 
       xlab="Eruptions (min)", ylab="Density", main=paste("Chain", CH), 
       breaks=seq(1.4, 5.6, by=0.3), ylim=ylimE)
  lines(MPDensModel0[[CH]]$x$x1, MPDensModel0[[CH]]$dens[[1]], col="darkblue", lwd=2)
  hist(Faithful$waiting, prob=TRUE, col="sandybrown", 
       xlab="Waiting (min)", ylab="Density", main="", ylim=ylimW)
  lines(MPDensModel0[[CH]]$x$x2, MPDensModel0[[CH]]$dens[[2]], col="darkblue", lwd=2)
}  
dev.off()


###################################################
### code chunk number 22: Model0: Contour plot of the joint predictive density with the scatterplot
###################################################
postscript(paste(FIGDIR, "figFaithful05.ps", sep=""), width=7, height=10, 
           horizontal=FALSE)
par(mfrow=c(2, 1), bty="n")
for (CH in 1:2){
  plot(Faithful, col="red", xlab="Eruptions (min)", ylab="Waiting (min)", 
       main=paste("Chain", CH))
  contour(JPDensModel0[[CH]]$x$x1, JPDensModel0[[CH]]$x$x2, 
          JPDensModel0[[CH]]$dens[["1-2"]], 
          col="darkblue", add=TRUE)
}  
dev.off()


###################################################
### code chunk number 23: Model0: Image plot of the joint predictive density with the scatterplot
###################################################
postscript(paste(FIGDIR, "figFaithful06.ps", sep=""), width=7, height=10, 
           horizontal=FALSE)
par(mfrow=c(2, 1), bty="n")
for (CH in 1:2){
  plot(Faithful, col="darkblue", xlab="Eruptions (min)", ylab="Waiting (min)", 
       main=paste("Chain", CH))
  image(JPDensModel0[[CH]]$x$x1, JPDensModel0[[CH]]$x$x2, 
        JPDensModel0[[CH]]$dens[["1-2"]], add=TRUE, 
        col=rev(heat_hcl(33, c=c(80, 30), l=c(30, 90), power=c(1/5, 1.3))))  
  points(Faithful, col="darkblue")
}  
dev.off()


###################################################
### code chunk number 24: RJMCMC: Choose chain
###################################################
CH <- 1


###################################################
### code chunk number 25: Model0: Create mcmc objects
###################################################
if (RUN.ALLOUT){
  start <- Model0[[CH]]$nMCMC["burn"] + 1
  end <- Model0[[CH]]$nMCMC["burn"] + Model0[[CH]]$nMCMC["keep"]
  chgammaInv <- mcmc(Model0[[CH]]$gammaInv, start=start, end=end)
  chmixture <- mcmc(Model0[[CH]]$mixture, start=start, end=end)
  chdeviance <- mcmc(Model0[[CH]]$deviance, start=start, end=end)
}  


###################################################
### code chunk number 26: Choose iters to draw traceplots
###################################################
if (RUN.ALLOUT){
  tstart <- 495001
  tend <- 500000
  titers <- tstart:tend
  chgammaInv2 <- mcmc(Model0[[CH]]$gammaInv[titers,], start=tstart, end=tend)
  chmixture2 <- mcmc(Model0[[CH]]$mixture[titers,], start=tstart, end=tend)
  chdeviance2 <- mcmc(Model0[[CH]]$deviance[titers,], start=tstart, end=tend)
}


###################################################
### code chunk number 27: Model0: Traceplots
###################################################
if (RUN.ALLOUT){
  lwd <- 0.5
  postscript(paste(FIGKEEPDIR, "figFaithful07.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(mfrow=c(2, 3), bty="n")
  traceplot(chmixture2[, "y.Mean.1"], smooth=FALSE, col="darkblue", lwd=lwd, 
            main="E(eruptions)")
  traceplot(chmixture2[, "y.Mean.2"], smooth=FALSE, col="darkblue", lwd=lwd, 
            main="E(waitings)")
  traceplot(chgammaInv2[, "gammaInv1"], smooth=FALSE, col="brown", lwd=lwd, 
            main="gamma^{-1}")
  traceplot(chmixture2[, "y.SD.1"], smooth=FALSE, col="darkgreen", lwd=lwd, 
            main="sd(eruptions)")
  traceplot(chmixture2[, "y.SD.2"], smooth=FALSE, col="darkgreen", lwd=lwd, 
            main="sd(waitings)")
  traceplot(chmixture2[, "y.Corr.2.1"], smooth=FALSE, col="red", lwd=lwd, 
            main="corr(eruptions, waitings)")
  dev.off()
}


###################################################
### code chunk number 28: Model0: Traceplots
###################################################
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figFaithful08.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(mfrow=c(2, 2), bty="n")
  traceplot(chdeviance2[, "LogL0"], smooth=FALSE, col="red", lwd=lwd, 
            main="Log(L0)")
  traceplot(chdeviance2[, "LogL1"], smooth=FALSE, col="red", lwd=lwd, 
            main="Log(L1)")
  traceplot(chdeviance2[, "dev.complete"], smooth=FALSE, col="red", lwd=lwd, 
            main="D(complete)")
  traceplot(chdeviance2[, "dev.observed"], smooth=FALSE, col="red", lwd=lwd, 
            main="D(observed)")
  dev.off()
}  


###################################################
### code chunk number 29: Model0: Density plots
###################################################
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figFaithful09.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(mfrow=c(2, 3), bty="n")
  densplot(chmixture[, "y.Mean.1"], show.obs=FALSE, col="darkblue", 
           main="E(eruptions)")
  densplot(chmixture[, "y.Mean.2"], show.obs=FALSE, col="darkblue", 
           main="E(waitings)")
  densplot(chgammaInv[, "gammaInv1"], show.obs=FALSE, col="brown", 
           main="gamma^{-1}")
  densplot(chmixture[, "y.SD.1"], show.obs=FALSE, col="darkgreen", 
           main="sd(eruptions)")
  densplot(chmixture[, "y.SD.2"], show.obs=FALSE, col="darkgreen", 
           main="sd(waitings)")
  densplot(chmixture[, "y.Corr.2.1"], show.obs=FALSE, col="red", 
           main="corr(eruptions, waitings)")
  dev.off()
} 


###################################################
### code chunk number 30: Model0: Autocorrelation plots
###################################################
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figFaithful10.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(mfrow=c(2, 3), bty="n")
  autocorr.plot(chmixture[, "y.Mean.1"], auto.layout=FALSE, 
                ask=FALSE, col="darkblue", main="E(eruptions)")
  autocorr.plot(chmixture[, "y.Mean.2"], auto.layout=FALSE, 
                ask=FALSE, col="darkblue", main="E(waitings)")
  autocorr.plot(chgammaInv[, "gammaInv1"], auto.layout=FALSE, 
                ask=FALSE, col="brown", main="gamma^{-1}")
  autocorr.plot(chmixture[, "y.SD.1"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="sd(eruptions)")
  autocorr.plot(chmixture[, "y.SD.2"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="sd(waitings)")
  autocorr.plot(chmixture[, "y.Corr.2.1"], auto.layout=FALSE, 
                ask=FALSE, col="red", main="corr(eruptions, waitings)")
  dev.off()
}


###################################################
### code chunk number 31: Model0: Scatterplots and histograms of mixture elements
###################################################
if (RUN.ALLOUT){
  PCH <- 1;  CEX <- 0.5
  set.seed(770328)
  SELECT <- sample(end-start+1, size=500, replace=FALSE)

  MuDens1 <- MuDens2 <- list()
  for (j in 1:3){
    MuDens1[[j]] <- density(Model0[[CH]]$mu[,(j-1)*2+1])
    MuDens2[[j]] <- density(Model0[[CH]]$mu[,j*2])
  }  
  LTY <- c(1, 2, 4)  
  COL <- c("blue", "red", "darkgreen")
  
  XLIM <- c(-2, 2);  YLIM <- c(-2, 2);  
  XLAB <- "Mean, margin 1";  YLAB <- "Mean, margin 2"  
  postscript(paste(FIGKEEPDIR, "figFaithful16.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(bty="n")
  layout(matrix(c(1,0, 1,4, 2,4, 2,5, 3,5, 3,6), ncol=2, byrow=TRUE))  
  plot(Model0[[CH]]$mu[SELECT, "mu.1.1"], Model0[[CH]]$mu[SELECT, "mu.1.2"], 
       col="blue", pch=PCH, cex=CEX, xlim=XLIM, ylim=YLIM, 
       xlab=XLAB, ylab=YLAB, main=expression(mu[1]))
  plot(Model0[[CH]]$mu[SELECT, "mu.2.1"], Model0[[CH]]$mu[SELECT, "mu.2.2"], 
       col="blue", pch=PCH, cex=CEX, xlim=XLIM, ylim=YLIM, 
       xlab=XLAB, ylab=YLAB, main=expression(mu[2]))
  plot(Model0[[CH]]$mu[SELECT, "mu.3.1"], Model0[[CH]]$mu[SELECT, "mu.3.2"], 
       col="blue", pch=PCH, cex=CEX, xlim=XLIM, ylim=YLIM, 
       xlab=XLAB, ylab=YLAB, main=expression(mu[3]))
  #
  plot(MuDens1[[1]]$x, MuDens1[[1]]$y, type="l", lty=LTY[1], col=COL[1], 
       xlab="Mean, margin 1", ylab="Posterior density", xlim=XLIM, 
       main="Margin 1")
  for (j in 2:3) lines(MuDens1[[j]]$x, MuDens1[[j]]$y, lty=LTY[j], col=COL[j])  
  #
  plot(MuDens2[[1]]$x, MuDens2[[1]]$y, type="l", lty=LTY[1], col=COL[1], 
       xlab="Mean, margin 2", ylab="Posterior density", xlim=YLIM, 
       main="Margin 2")
  for (j in 2:3) lines(MuDens2[[j]]$x, MuDens2[[j]]$y, lty=LTY[j], col=COL[j])
  #
  plot(0:100, 0:100, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
  legend(10, 99, paste("Comp.", 1:3), lty=LTY, col=COL, bty="n", y.intersp=1)    
  dev.off()
  
  XLIM <- c(0, 1)
  YLIM <- c(0, 1.5)
  XLAB <- "Variance, margin 1"
  YLAB <- "Variance, margin 2"   
  postscript(paste(FIGKEEPDIR, "figFaithful17.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(bty="n")
  layout(matrix(c(0,1,1,0, 2,2,3,3), nrow=2, byrow=TRUE))
  plot(Model0[[CH]]$Sigma[SELECT, "Sigma1.1.1"], 
       Model0[[CH]]$Sigma[SELECT, "Sigma1.2.2"], 
       col="red", pch=PCH, cex=CEX, xlim=XLIM, ylim=YLIM, 
       xlab=XLAB, ylab=YLAB, main=expression(Sigma[1]))
  plot(Model0[[CH]]$Sigma[SELECT, "Sigma2.1.1"], 
       Model0[[CH]]$Sigma[SELECT, "Sigma2.2.2"], 
       col="red", pch=PCH, cex=CEX, xlim=XLIM, ylim=YLIM, 
       xlab=XLAB, ylab=YLAB, main=expression(Sigma[2]))
  plot(Model0[[CH]]$Sigma[SELECT, "Sigma3.1.1"], 
       Model0[[CH]]$Sigma[SELECT, "Sigma3.2.2"], 
       col="red", pch=PCH, cex=CEX, xlim=XLIM, ylim=YLIM, 
       xlab=XLAB, ylab=YLAB, main=expression(Sigma[3]))
  dev.off()
  
  XLIM <- c(0, 1)
  YLIM <- c(0, 4.5)
  XLAB <- "Weight"
  YLAB <- "Density"
  postscript(paste(FIGKEEPDIR, "figFaithful18.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(bty="n")
  layout(matrix(c(0,1,1,0, 2,2,3,3), nrow=2, byrow=TRUE))
  hist(Model0[[CH]]$w[, "w1"], prob=TRUE, col="sandybrown", 
       xlim=XLIM, ylim=YLIM, xlab=XLAB, ylab=YLAB, main=expression(w[1]))  
  hist(Model0[[CH]]$w[, "w2"], prob=TRUE, col="sandybrown", 
       xlim=XLIM, ylim=YLIM, xlab=XLAB, ylab=YLAB, main=expression(w[2]))  
  hist(Model0[[CH]]$w[, "w3"], prob=TRUE, col="sandybrown", 
       xlim=XLIM, ylim=YLIM, xlab=XLAB, ylab=YLAB, main=expression(w[3]))  
  dev.off()
}  


###################################################
### code chunk number 32: Fixed K MCMC
###################################################
if (RUN.TIMECONSUMING.CODE){
  Seed <- c(777988621, 777988621, 777988621, 777988621, 780830,
            780830, 777988621, 777988621, 777988621, 777988621)
#
  Keep <- c("iter", "nMCMC", "dim", "prior", "init", "RJMCMC", 
            "scale", "state", "freqK", "propK", "DIC", "moves", 
            "pm.y", "pm.z", "pm.indDev", "pred.dens", "summ.y.Mean", 
            "summ.y.SDCorr", "summ.z.Mean", "summ.z.SDCorr")
#
  ModelK <- list()
  MPDensModelK <- list()
  JPDensModelK <- list()
  for (k in 1:10){
    set.seed(Seed[k])    
    cat(paste("K = ", k, "\n-------------------------------\n", sep=""))  
    PriorNow <- Prior0
    PriorNow$Kmax <- k
    ModelK[[k]] <- NMixMCMC(y0=Faithful, prior=PriorNow, nMCMC=nMCMC, PED=TRUE)
#
    cat(paste("\nComputation of marginal pred. densities started on ", date(), 
              "\n", sep=""))  
    MPDensModelK[[k]] <- list()
    MPDensModelK[[k]][[1]] <- NMixPredDensMarg(ModelK[[k]][[1]], grid=ygrid)
    MPDensModelK[[k]][[2]] <- NMixPredDensMarg(ModelK[[k]][[2]], grid=ygrid)
    cat(paste("Computation of marginal pred. densities finished on ", date(), 
              "\n\n", sep=""))  
#
    cat(paste("Computation of joint pred. densities started on ", date(), 
              "\n", sep=""))    
    JPDensModelK[[k]] <- list()
    JPDensModelK[[k]][[1]] <- NMixPredDensJoint2(ModelK[[k]][[1]], grid=ygrid)
    JPDensModelK[[k]][[2]] <- NMixPredDensJoint2(ModelK[[k]][[2]], grid=ygrid)
    cat(paste("Computation of joint pred. densities finished on ", date(), 
              "\n\n\n", sep=""))      
#
    ModelK[[k]][[1]] <- ModelK[[k]][[1]][Keep]
    ModelK[[k]][[2]] <- ModelK[[k]][[2]][Keep]  
    class(ModelK[[k]][[1]]) <- class(ModelK[[k]][[2]]) <- "NMixMCMC"
  }
}  


###################################################
### code chunk number 33: Fixed K: PED and DIC
###################################################
PED <- ModelK[[1]]$PED
DIC <- list(Chain1=ModelK[[1]][[1]]$DIC, Chain2=ModelK[[1]][[2]]$DIC)
for (k in 2:length(ModelK)){
  PED <- rbind(PED, ModelK[[k]]$PED)
  DIC[[1]] <- rbind(DIC[[1]], ModelK[[k]][[1]]$DIC)
  DIC[[2]] <- rbind(DIC[[2]], ModelK[[k]][[2]]$DIC)  
}
rownames(PED) <- rownames(DIC[[1]]) <- rownames(DIC[[2]]) <- paste("K=", 1:length(ModelK), sep="")


###################################################
### code chunk number 34: Fixed K: print PED
###################################################
print(PED)


###################################################
### code chunk number 35: Fixed K: print DIC
###################################################
print(DIC)


###################################################
### code chunk number 36: Fixed K MCMC: Eruptions - plot of the marginal predictive density
###################################################
CH <- 1
postscript(paste(FIGDIR, "figFaithful11.ps", sep=""), width=6, height=6, 
           horizontal=FALSE)
par(mfrow=c(1, 1), bty="n")
hist(Faithful$eruptions, prob=TRUE, col="grey90", 
     xlab="Eruptions (min)", ylab="Density", main="", 
     breaks=seq(1.4, 5.6, by=0.3))
for (k in 1:10){
  lines(MPDensModelK[[k]][[CH]]$x$x1, MPDensModelK[[k]][[CH]]$dens[[1]], 
        col="red")
}
dev.off()


###################################################
### code chunk number 37: Fixed K MCMC: Eruptions - plot of the marginal predictive density
###################################################
postscript(paste(FIGDIR, "figFaithful12.ps", sep=""), width=7, height=10, 
           horizontal=FALSE)
par(mar=c(3, 2, 2, 1)+0.1)
par(mfrow=c(5, 2), bty="n")
for (k in 1:10){
  hist(Faithful$eruptions, prob=TRUE, col="lightblue", 
       xlab="", ylab="", main=paste("K = ", k, sep=""), 
       breaks=seq(1.4, 5.6, by=0.3))
  lines(MPDensModelK[[k]][[CH]]$x$x1, MPDensModelK[[k]][[CH]]$dens[[1]], 
        col="red", lwd=2)
}
dev.off()


###################################################
### code chunk number 38: Fixed K MCMC: Waiting - plot of the marginal predictive density
###################################################
CH <- 1
postscript(paste(FIGDIR, "figFaithful13.ps", sep=""), width=6, height=6, 
           horizontal=FALSE)
par(mfrow=c(1, 1), bty="n")
hist(Faithful$waiting, prob=TRUE, col="grey90", 
     xlab="Waiting (min)", ylab="Density", main="")
for (k in 1:10){
  lines(MPDensModelK[[k]][[CH]]$x$x2, MPDensModelK[[k]][[CH]]$dens[[2]], 
        col="red")
}
dev.off()


###################################################
### code chunk number 39: Fixed K MCMC: Waiting - plot of the marginal predictive density
###################################################
postscript(paste(FIGDIR, "figFaithful14.ps", sep=""), width=7, height=10, 
           horizontal=FALSE)
par(mar=c(3, 2, 2, 1)+0.1)
par(mfrow=c(5, 2), bty="n")
for (k in 1:10){
  hist(Faithful$waiting, prob=TRUE, col="lightblue", 
       xlab="", ylab="", main=paste("K = ", k, sep=""))
  lines(MPDensModelK[[k]][[CH]]$x$x2, MPDensModelK[[k]][[CH]]$dens[[2]], 
        col="red", lwd=2)
}  
dev.off()


###################################################
### code chunk number 40: Fixed K MCMC: Waiting - plot of the joint predictive density
###################################################
CH <- 1
postscript(paste(FIGDIR, "figFaithful15.ps", sep=""), width=7, height=10, 
           horizontal=FALSE)
par(mar=c(3, 2, 2, 1)+0.1)
par(mfrow=c(5, 2), bty="n")
for (k in 1:10){
  plot(Faithful, col="red", xlab="", ylab="", 
       main=paste("K = ", k, sep=""))
  contour(JPDensModelK[[k]][[CH]]$x$x1, JPDensModelK[[k]][[CH]]$x$x2, 
          JPDensModelK[[k]][[CH]]$dens[["1-2"]], 
          col="darkblue", add=TRUE)
}  
dev.off()


###################################################
### code chunk number 41: Fixed K MCMC: Waiting - plot of the joint predictive density
###################################################
CH <- 1
postscript(paste(FIGDIR, "figFaithful19.ps", sep=""), width=7, height=10, 
           horizontal=FALSE)
par(mar=c(3, 2, 2, 1)+0.1)
par(mfrow=c(5, 2), bty="n")
for (k in 1:10){
  plot(Faithful, col="darkblue", xlab="", ylab="", 
       main=paste("K = ", k, sep=""))
  image(JPDensModelK[[k]][[CH]]$x$x1, JPDensModelK[[k]][[CH]]$x$x2, 
        JPDensModelK[[k]][[CH]]$dens[["1-2"]], add=TRUE,
        col=rev(heat_hcl(33, c=c(80, 30), l=c(30, 90), power=c(1/5, 1.3))))  
#  points(Faithful, col="darkblue")
}  
dev.off()


###################################################
### code chunk number 42: Save results
###################################################
if (RUN.TIMECONSUMING.CODE){
  save(list="Model0", 
       file=paste(RESULTDIR, "/Faithful-Model0", Kshow, ".RData", sep=""))        
  save(list=c("ModelK", "MPDensModelK", "JPDensModelK"), 
       file=paste(RESULT2DIR, "/Faithful-Result.RData", sep=""))    
}  


