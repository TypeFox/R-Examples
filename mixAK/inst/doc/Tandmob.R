###
###  Full pdf document describing the code included here is available at
###  http://msekce.karlin.mff.cuni.cz/~komarek/software/mixAK/Tandmob.pdf
###
### ==============================================================================

###################################################
### chunk number 1: What should be created in this Sweave run
###################################################
#line 89 "Tandmob.Rnw"
RUN.TIMECONSUMING.CODE <- FALSE
RUN.ALLOUT <- FALSE


###################################################
### chunk number 2: Directory to store figures
###################################################
#line 97 "Tandmob.Rnw"
FIGDIR <- "./figures/"
FIGKEEPDIR <- "./figuresKeep/"


###################################################
### chunk number 3: Directory with results computed in past
###################################################
#line 108 "Tandmob.Rnw"
RESULTDIR <- "/home/komarek/RESULT_OBJ/mixAK-Tandmob-S090426/"  ### must be changed by the user
RESULT2DIR <- "./RESULT_OBJ/"      ### must be changed by the user


###################################################
### chunk number 4: Options
###################################################
#line 114 "Tandmob.Rnw"
options(width=80)


###################################################
### chunk number 5: Check for existence of directories to store results
###################################################
#line 119 "Tandmob.Rnw"
if (!file.exists(RESULTDIR)){
  stop(paste("Directory ", RESULTDIR, " does not exist.\nYou have to create it or change the value of the variable RESULTDIR.\n"))
}  
if (!file.exists(RESULT2DIR)){
  stop(paste("Directory ", RESULT2DIR, " does not exist.\nYou have to create it or change the value of the variable RESULT2DIR.\n"))
}  


###################################################
### chunk number 6: Load results computed in past
###################################################
#line 131 "Tandmob.Rnw"
Kshow <- 2

if (("Tandmob-Result.RData" %in% dir(RESULT2DIR)) & (paste("Tandmob-PDensBiModelK0", Kshow,".RData", sep="") %in% dir(RESULT2DIR))){
  load(paste(RESULT2DIR, "Tandmob-Result.RData", sep=""))                    
     ## contains ModelK (without chains), PDensUniModelK
  load(paste(RESULT2DIR, "Tandmob-PDensBiModelK0", Kshow,".RData", sep=""))  
     ## contains PDensBiModelK[[Kshow]]  
  Model0 <- ModelK[[Kshow]]
  PDensUniModel0 <- PDensUniModelK[[Kshow]]
  PDensBiModel0 <- PDensBiModelK[[Kshow]]  
}else{
  if (!RUN.TIMECONSUMING.CODE){    
    stop(paste("Directory ", RESULT2DIR, " does not contain necessary files.\nSet RUN.TIMECONSUMING.CODE to TRUE.", sep=""))
  }        
}  

if (RUN.ALLOUT){
  if (paste("Tandmob-Model0", Kshow, ".RData", sep="") %in% dir(RESULTDIR)){
    load(paste(RESULTDIR, "Tandmob-Model0", Kshow, ".RData", sep=""))          
       ## contains Model0=ModelK[[Kshow]] (chains included)  
  }else{      
    stop(paste("Directory ", RESULTDIR, " does not contain necessary files.\nSet RUN.TIMECONSUMING.CODE to TRUE.", sep=""))
  }  
} 


###################################################
### chunk number 7: Load needed packages
###################################################
#line 161 "Tandmob.Rnw"
library("mixAK")
library("coda")
library("colorspace")


###################################################
### chunk number 8: Read the data and compute a brief summary
###################################################
#line 173 "Tandmob.Rnw"
data("Tandmob", package="mixAK")
data("TandmobEmer", package="mixAK")


###################################################
### chunk number 9: Table showing the number and proportion of each type of censoring for each tooth
###################################################
#line 178 "Tandmob.Rnw"
NUM.PROP.CENS <- function(data, tanden)
{
  TABLE <- matrix(NA, nrow=length(tanden), ncol=3)
  rownames(TABLE) <- paste(tanden)
  colnames(TABLE) <- c("Left", "Interval", "Right")
  for (tt in 1:length(tanden)){
    ebeg <- get(data)[,paste("EBEG.", tanden[tt], sep="")]
    eend <- get(data)[,paste("EEND.", tanden[tt], sep="")]  
    TABLE[tt, "Left"]     <- sum(is.na(ebeg) & !is.na(eend))  
    TABLE[tt, "Interval"] <- sum(!is.na(ebeg) & !is.na(eend))
    TABLE[tt, "Right"]    <- sum(!is.na(ebeg) & is.na(eend))      
    rm(list=c("ebeg", "eend"))  
  }
  PROP.TABLE <- t(round(apply(TABLE, 1, prop.table)*100, 3))
  RET <- list(Table=TABLE, Prop.Table=PROP.TABLE)
  return(RET)
}

tanden <- rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4)
TAB01 <- NUM.PROP.CENS(data="Tandmob", tanden=tanden)


###################################################
### chunk number 10: Print table TAB01
###################################################
#line 202 "Tandmob.Rnw"
print(TAB01)


###################################################
### chunk number 11: Preparation of the data using the original data.frame Tandmob
###################################################
#line 221 "Tandmob.Rnw"
Emerg.min <- 5

Tooth <- 10 + 1:6
y0 <- Tandmob[, paste("EBEG.", Tooth, sep="")]
y1 <- Tandmob[, paste("EEND.", Tooth, sep="")]
y0[is.na(y0)] <- Emerg.min     
   ### Left-censored changed into interval-censored with the lower limit=5
   ###
censor <- matrix(3, nrow=nrow(y0), ncol=length(Tooth))
censor[is.na(y1)] <- 0
   ###
colnames(y0) <- colnames(y1) <- colnames(censor) <- Tooth
rownames(y0) <- rownames(y1) <- rownames(censor) <- Tandmob$IDNR


###################################################
### chunk number 12: Preparation of the data using the data.frame TandmobEmer
###################################################
#line 239 "Tandmob.Rnw"
Tooth <- 10 + 1:6
y0 <- TandmobEmer[, paste("EBEG.", Tooth, sep="")]
y1 <- TandmobEmer[, paste("EEND.", Tooth, sep="")]
censor <- TandmobEmer[, paste("CENSOR.", Tooth, sep="")]
colnames(y0) <- colnames(y1) <- colnames(censor) <- Tooth
rownames(y0) <- rownames(y1) <- rownames(censor) <- TandmobEmer$IDNR


###################################################
### chunk number 13: Print part of prepared data
###################################################
#line 251 "Tandmob.Rnw"
print(y0[1:5,])
print(y1[1:5,])


###################################################
### chunk number 14: Print censor indicators
###################################################
#line 258 "Tandmob.Rnw"
print(censor[1:5,])


###################################################
### chunk number 15: Length of the MCMC
###################################################
#line 265 "Tandmob.Rnw"
nMCMC <- c(burn=10000, keep=20000, thin=10, info=1000)


###################################################
### chunk number 16: Grid of values for the predictive density
###################################################
#line 272 "Tandmob.Rnw"
lygridBi <- 50
lygridUni <- 100
ymin <- c(4, 5, 6, 6, 6, 4)
ymax <- c(10, 12, 17, 16, 17, 9)
ygridBi <- list(seq(ymin[1], ymax[1], length=lygridBi),
                seq(ymin[2], ymax[2], length=lygridBi),
                seq(ymin[3], ymax[3], length=lygridBi),
                seq(ymin[4], ymax[4], length=lygridBi),
                seq(ymin[5], ymax[5], length=lygridBi),
                seq(ymin[6], ymax[6], length=lygridBi))
ygridUni <- list(seq(ymin[1], ymax[1], length=lygridUni),
                 seq(ymin[2], ymax[2], length=lygridUni),
                 seq(ymin[3], ymax[3], length=lygridUni),
                 seq(ymin[4], ymax[4], length=lygridUni),
                 seq(ymin[5], ymax[5], length=lygridUni),
                 seq(ymin[6], ymax[6], length=lygridUni))


###################################################
### chunk number 17: TESTING: Prior distribution and model with two components and posterior predictive density
###################################################
#line 299 "Tandmob.Rnw"
set.seed(770328)
Prior0 <- list(priorK="fixed", Kmax=2)
TEST.Model0 <- NMixMCMC(y0=y0[1:100,], y1=y1[1:100,], censor=censor[1:100,], prior=Prior0, 
                        nMCMC=c(burn=100, keep=100, thin=2, info=50), scale=list(shift=0, scale=1), PED=TRUE)
TEST.Model0 <- NMixMCMC(y0=y0, y1=y1, censor=censor, prior=Prior0, 
                        nMCMC=c(burn=10, keep=10, thin=2, info=10), scale=list(shift=0, scale=1), PED=TRUE)
TEST.PDensUniModel0 <- NMixPredDensMarg(TEST.Model0[[1]], grid=ygridUni)
TEST.PDensBiModel0 <- NMixPredDensJoint2(TEST.Model0[[1]], grid=ygridBi)


###################################################
### chunk number 18: Prior distribution
###################################################
#line 312 "Tandmob.Rnw"
Prior0 <- list(priorK="fixed", Kmax=2)


###################################################
### chunk number 19: Model with two components
###################################################
#line 318 "Tandmob.Rnw"
if (RUN.TIMECONSUMING.CODE){
  set.seed(770328)  
  Model0 <- NMixMCMC(y0=y0, y1=y1, censor=censor, prior=Prior0, 
                     nMCMC=nMCMC, scale=list(shift=0, scale=1), PED=TRUE)
}


###################################################
### chunk number 20: The prior distribution for the function NMixMCMC was the same as with eval=FALSE
###################################################
## #line 349 "Tandmob.Rnw"
## Prior0 <- list(priorK="fixed", Kmax=2,
##             delta=1,
##             priormuQ="independentC", 
##             xi=c(8.425, 9.508172, 9.769336, 9.75195, 9.789845, 7.675), 
##             D=diag(c(31.9225, 59.54195, 67.87572, 67.30397, 68.55327, 18.0625)),
##             zeta=7, 
##             g=0.2, 
##             h=c(0.3132587, 0.1679488, 0.1473281, 0.1488796, 0.1458720, 0.5536332))


###################################################
### chunk number 21: Model0: Basic posterior summary
###################################################
#line 396 "Tandmob.Rnw"
print(Model0)


###################################################
### chunk number 22: Model0: Marginal predictive densities
###################################################
#line 403 "Tandmob.Rnw"
if (RUN.TIMECONSUMING.CODE){
  PDensUniModel0 <- list()
  PDensUniModel0[[1]] <- NMixPredDensMarg(Model0[[1]], grid=ygridUni)  
  PDensUniModel0[[2]] <- NMixPredDensMarg(Model0[[2]], grid=ygridUni)  
}


###################################################
### chunk number 23: Model0: Plot of the marginal predictive density
###################################################
#line 413 "Tandmob.Rnw"
postscript(paste(FIGDIR, "figTandmob01.ps", sep=""), width=10, height=7, 
           horizontal=FALSE)
plot(PDensUniModel0[[1]])
dev.off()


###################################################
### chunk number 24: Model0: Joint predictive density
###################################################
#line 433 "Tandmob.Rnw"
if (RUN.TIMECONSUMING.CODE){
  PDensBiModel0 <- list()
  PDensBiModel0[[1]] <- NMixPredDensJoint2(Model0[[1]], grid=ygridBi)  
  PDensBiModel0[[2]] <- NMixPredDensJoint2(Model0[[2]], grid=ygridBi)  
}


###################################################
### chunk number 25: Model0: Plot of the joint predictive density
###################################################
#line 443 "Tandmob.Rnw"
postscript(paste(FIGDIR, "figTandmob02.ps", sep=""), width=20, height=14, 
           horizontal=FALSE)
plot(PDensBiModel0[[1]])
dev.off()


###################################################
### chunk number 26: Model0: Nicer plot of the marginal predictive density
###################################################
#line 466 "Tandmob.Rnw"
postscript(paste(FIGDIR, "figTandmob03.ps", sep=""), width=7, height=5, 
           horizontal=FALSE)
xlim <- c(4, 17)
ylim <- c(0, 0.8)
par(mfrow=c(2, 3), bty="n", mar=c(4, 4, 4, 0)+0.1)
for (tt in 1:6){
  plot(PDensUniModel0[[1]]$x[[tt]], PDensUniModel0[[1]]$dens[[tt]], 
       type="l", col="red", xlim=xlim, ylim=ylim, 
       xlab="Age (years)", ylab="Density of emergence", 
       main=paste("Tooth ", 10+tt, sep=""))
  lines(PDensUniModel0[[2]]$x[[tt]], PDensUniModel0[[2]]$dens[[tt]], col="blue")
}  
dev.off()


###################################################
### chunk number 27: Model0: Contour plot of the joint bivariate predictive densities with the scatterplot
###################################################
#line 485 "Tandmob.Rnw"
CH <- 1
postscript(paste(FIGDIR, "figTandmob04.ps", sep=""), width=10, height=7, 
           horizontal=FALSE)
par(mfrow=c(3, 5), bty="n", mar=c(4, 4, 4, 0)+0.1)
for (tt1 in 1:5){
  for (tt2 in (tt1+1):6){
    contour(PDensBiModel0[[CH]]$x[[tt1]], PDensBiModel0[[CH]]$x[[tt2]], 
            PDensBiModel0[[CH]]$dens[[paste(tt1, "-", tt2, sep="")]],
            col="red", xlab="Age (years)", ylab="Age (years)", 
            main=paste("Teeth ", 10+tt1, "-", 10+tt2, sep=""))    
  }  
}
dev.off()


###################################################
### chunk number 28: Model0: Image plot of the joint bivariate predictive densities with the scatterplot
###################################################
#line 522 "Tandmob.Rnw"
CH <- 1
postscript(paste(FIGDIR, "figTandmob05.ps", sep=""), width=10, height=7, 
           horizontal=FALSE)
par(mfrow=c(3, 5), bty="n", mar=c(4, 4, 4, 0)+0.1)
for (tt1 in 1:5){
  for (tt2 in (tt1+1):6){
    image(PDensBiModel0[[CH]]$x[[tt1]], PDensBiModel0[[CH]]$x[[tt2]], 
          PDensBiModel0[[CH]]$dens[[paste(tt1, "-", tt2, sep="")]],
          col=rev(heat_hcl(33, c=c(80, 30), l=c(30, 90), power=c(1/5, 1.3))),
          xlab="Age (years)", ylab="Age (years)", 
          main=paste("Teeth ", 10+tt1, "-", 10+tt2, sep=""))    
  }  
}
dev.off()


###################################################
### chunk number 29: RJMCMC: Choose chain
###################################################
#line 554 "Tandmob.Rnw"
CH <- 1


###################################################
### chunk number 30: Model0: Create mcmc objects
###################################################
#line 559 "Tandmob.Rnw"
if (RUN.ALLOUT){
  start <- Model0[[CH]]$nMCMC["burn"] + 1
  end <- Model0[[CH]]$nMCMC["burn"] + Model0[[CH]]$nMCMC["keep"]
  chgammaInv <- mcmc(Model0[[CH]]$gammaInv, start=start, end=end)
  chmixture <- mcmc(Model0[[CH]]$mixture, start=start, end=end)
  chdeviance <- mcmc(Model0[[CH]]$deviance, start=start, end=end)
}  


###################################################
### chunk number 31: Model0: Traceplots
###################################################
#line 573 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob06.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)  
  lwd <- 0.5
  par(mfrow=c(2, 3), bty="n")
  traceplot(chmixture[, "y.Mean.1"], smooth=FALSE, col="darkblue", 
            lwd=lwd, main="E(Tooth 11)")
  traceplot(chmixture[, "y.Mean.2"], smooth=FALSE, col="darkblue", 
            lwd=lwd, main="E(Tooth 12)")
  traceplot(chmixture[, "y.Mean.3"], smooth=FALSE, col="darkblue", 
            lwd=lwd, main="E(Tooth 13)")
  traceplot(chmixture[, "y.Mean.4"], smooth=FALSE, col="darkblue", 
            lwd=lwd, main="E(Tooth 14)")
  traceplot(chmixture[, "y.Mean.5"], smooth=FALSE, col="darkblue", 
            lwd=lwd, main="E(Tooth 15)")
  traceplot(chmixture[, "y.Mean.6"], smooth=FALSE, col="darkblue", 
            lwd=lwd, main="E(Tooth 16)")
  dev.off()
}  


###################################################
### chunk number 32: Model0: Traceplots
###################################################
#line 607 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob07.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(mfrow=c(2, 3), bty="n")
  traceplot(chmixture[, "y.SD.1"], smooth=FALSE, col="darkgreen", 
            lwd=lwd, main="SD(Tooth 11)")
  traceplot(chmixture[, "y.SD.2"], smooth=FALSE, col="darkgreen", 
            lwd=lwd, main="SD(Tooth 12)")
  traceplot(chmixture[, "y.SD.3"], smooth=FALSE, col="darkgreen", 
            lwd=lwd, main="SD(Tooth 13)")
  traceplot(chmixture[, "y.SD.4"], smooth=FALSE, col="darkgreen", 
            lwd=lwd, main="SD(Tooth 14)")
  traceplot(chmixture[, "y.SD.5"], smooth=FALSE, col="darkgreen", 
            lwd=lwd, main="SD(Tooth 15)")
  traceplot(chmixture[, "y.SD.6"], smooth=FALSE, col="darkgreen", 
            lwd=lwd, main="SD(Tooth 16)")
  dev.off()
}


###################################################
### chunk number 33: Model0: Traceplots
###################################################
#line 630 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob08.ps", sep=""), width=10, height=7, 
             horizontal=FALSE)
  par(bty="n", mar=c(1, 1, 4, 0)+0.1)
  layout(matrix(c(1,2,3,4,5, 0,6,7,8,9, 0,0,10,11,12, 0,0,0,13,14, 0,0,0,0,15), 
                nrow=5, byrow=TRUE))
  traceplot(chmixture[, "y.Corr.2.1"], smooth=FALSE, col="darkgreen", 
            main="Corr 2-1")
  traceplot(chmixture[, "y.Corr.3.1"], smooth=FALSE, col="darkgreen", 
            main="Corr 3-1")
  traceplot(chmixture[, "y.Corr.4.1"], smooth=FALSE, col="darkgreen", 
            main="Corr 4-1")
  traceplot(chmixture[, "y.Corr.5.1"], smooth=FALSE, col="darkgreen", 
            main="Corr 5-1")
  traceplot(chmixture[, "y.Corr.6.1"], smooth=FALSE, col="darkgreen", 
            main="Corr 6-1")
  traceplot(chmixture[, "y.Corr.3.2"], smooth=FALSE, col="darkgreen", 
            main="Corr 3-2")
  traceplot(chmixture[, "y.Corr.4.2"], smooth=FALSE, col="darkgreen", 
            main="Corr 4-2")
  traceplot(chmixture[, "y.Corr.5.2"], smooth=FALSE, col="darkgreen", 
            main="Corr 5-2")
  traceplot(chmixture[, "y.Corr.6.2"], smooth=FALSE, col="darkgreen", 
            main="Corr 6-2")
  traceplot(chmixture[, "y.Corr.4.3"], smooth=FALSE, col="darkgreen", 
            main="Corr 4-3")
  traceplot(chmixture[, "y.Corr.5.3"], smooth=FALSE, col="darkgreen", 
            main="Corr 5-3")
  traceplot(chmixture[, "y.Corr.6.3"], smooth=FALSE, col="darkgreen", 
            main="Corr 6-3")
  traceplot(chmixture[, "y.Corr.5.4"], smooth=FALSE, col="darkgreen", 
            main="Corr 5-4")
  traceplot(chmixture[, "y.Corr.6.4"], smooth=FALSE, col="darkgreen", 
            main="Corr 6-4")
  traceplot(chmixture[, "y.Corr.6.5"], smooth=FALSE, col="darkgreen", 
            main="Corr 6-5")
  dev.off()
}


###################################################
### chunk number 34: Model0: Traceplots
###################################################
#line 673 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob09.ps", sep=""), width=10, height=7, 
             horizontal=FALSE)
  par(bty="n")
  layout(matrix(c(1,1,2,2,3,3, 0,4,4,5,5,0), nrow=2, byrow=TRUE))
  traceplot(chgammaInv[, "gammaInv1"], smooth=FALSE, 
            col="brown", main="gamma^{-1}")
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
### chunk number 35: Model0: Posterior density estimates
###################################################
#line 719 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob10.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(mfrow=c(2, 3), bty="n")
  densplot(chmixture[, "y.Mean.1"], show.obs=FALSE, col="darkblue", 
           main="E(Tooth 11)")
  densplot(chmixture[, "y.Mean.2"], show.obs=FALSE, col="darkblue", 
           main="E(Tooth 12)")
  densplot(chmixture[, "y.Mean.3"], show.obs=FALSE, col="darkblue", 
           main="E(Tooth 13)")
  densplot(chmixture[, "y.Mean.4"], show.obs=FALSE, col="darkblue", 
           main="E(Tooth 14)")
  densplot(chmixture[, "y.Mean.5"], show.obs=FALSE, col="darkblue", 
           main="E(Tooth 15)")
  densplot(chmixture[, "y.Mean.6"], show.obs=FALSE, col="darkblue", 
           main="E(Tooth 16)")
  dev.off()
}


###################################################
### chunk number 36: Model0: Posterior density estimates
###################################################
#line 742 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob11.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(mfrow=c(2, 3), bty="n")
  densplot(chmixture[, "y.SD.1"], show.obs=FALSE, col="darkgreen",  
           main="SD(Tooth 11)")
  densplot(chmixture[, "y.SD.2"], show.obs=FALSE, col="darkgreen", 
           main="SD(Tooth 12)")
  densplot(chmixture[, "y.SD.3"], show.obs=FALSE, col="darkgreen", 
           main="SD(Tooth 13)")
  densplot(chmixture[, "y.SD.4"], show.obs=FALSE, col="darkgreen", 
           main="SD(Tooth 14)")
  densplot(chmixture[, "y.SD.5"], show.obs=FALSE, col="darkgreen", 
           main="SD(Tooth 15)")
  densplot(chmixture[, "y.SD.6"], show.obs=FALSE, col="darkgreen", 
           main="SD(Tooth 16)")
  dev.off()
}


###################################################
### chunk number 37: Model0: Posterior density estimates
###################################################
#line 781 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob12.ps", sep=""), width=10, height=7, 
             horizontal=FALSE)
  par(bty="n", mar=c(2, 1, 4, 0)+0.1)
  layout(matrix(c(1,2,3,4,5, 0,6,7,8,9, 0,0,10,11,12, 
                  0,0,0,13,14, 0,0,0,0,15), nrow=5, byrow=TRUE))
  densplot(chmixture[, "y.Corr.2.1"], show.obs=FALSE, col="darkgreen", 
           main="Corr 2-1")
  densplot(chmixture[, "y.Corr.3.1"], show.obs=FALSE, col="darkgreen", 
           main="Corr 3-1")
  densplot(chmixture[, "y.Corr.4.1"], show.obs=FALSE, col="darkgreen", 
           main="Corr 4-1")
  densplot(chmixture[, "y.Corr.5.1"], show.obs=FALSE, col="darkgreen", 
           main="Corr 5-1")
  densplot(chmixture[, "y.Corr.6.1"], show.obs=FALSE, col="darkgreen", 
           main="Corr 6-1")
  densplot(chmixture[, "y.Corr.3.2"], show.obs=FALSE, col="darkgreen", 
           main="Corr 3-2")
  densplot(chmixture[, "y.Corr.4.2"], show.obs=FALSE, col="darkgreen", 
           main="Corr 4-2")
  densplot(chmixture[, "y.Corr.5.2"], show.obs=FALSE, col="darkgreen", 
           main="Corr 5-2")
  densplot(chmixture[, "y.Corr.6.2"], show.obs=FALSE, col="darkgreen", 
           main="Corr 6-2")
  densplot(chmixture[, "y.Corr.4.3"], show.obs=FALSE, col="darkgreen", 
           main="Corr 4-3")
  densplot(chmixture[, "y.Corr.5.3"], show.obs=FALSE, col="darkgreen", 
           main="Corr 5-3")
  densplot(chmixture[, "y.Corr.6.3"], show.obs=FALSE, col="darkgreen", 
           main="Corr 6-3")
  densplot(chmixture[, "y.Corr.5.4"], show.obs=FALSE, col="darkgreen", 
           main="Corr 5-4")
  densplot(chmixture[, "y.Corr.6.4"], show.obs=FALSE, col="darkgreen", 
           main="Corr 6-4")
  densplot(chmixture[, "y.Corr.6.5"], show.obs=FALSE, col="darkgreen", 
           main="Corr 6-5")
  dev.off()
}


###################################################
### chunk number 38: Model0: Posterior density estimates
###################################################
#line 833 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob13.ps", sep=""), width=10, height=7, 
             horizontal=FALSE)
  par(bty="n")
  layout(matrix(c(1,1,2,2,3,3, 0,4,4,5,5,0), 
         nrow=2, byrow=TRUE))
  densplot(chgammaInv[, "gammaInv1"], show.obs=FALSE, 
           col="brown", main="gamma^{-1}")
  densplot(chdeviance[, "LogL0"], show.obs=FALSE, 
           col="red", main="Log(L0)")
  densplot(chdeviance[, "LogL1"], show.obs=FALSE, 
           col="red", main="Log(L1)")
  densplot(chdeviance[, "dev.complete"], show.obs=FALSE, 
           col="red", main="D(complete)")
  densplot(chdeviance[, "dev.observed"], show.obs=FALSE, 
           col="red", main="D(observed)")
  dev.off()
}


###################################################
### chunk number 39: Model0: Autocorrelation plots
###################################################
#line 866 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob14.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(mfrow=c(2, 3), bty="n")
  autocorr.plot(chmixture[, "y.Mean.1"], auto.layout=FALSE, ask=FALSE, 
                col="darkblue", main="E(Tooth 11)")
  autocorr.plot(chmixture[, "y.Mean.2"], auto.layout=FALSE, ask=FALSE, 
                col="darkblue", main="E(Tooth 12)")
  autocorr.plot(chmixture[, "y.Mean.3"], auto.layout=FALSE, ask=FALSE, 
                col="darkblue", main="E(Tooth 13)")
  autocorr.plot(chmixture[, "y.Mean.4"], auto.layout=FALSE, ask=FALSE, 
                col="darkblue", main="E(Tooth 14)")
  autocorr.plot(chmixture[, "y.Mean.5"], auto.layout=FALSE, ask=FALSE, 
                col="darkblue", main="E(Tooth 15)")
  autocorr.plot(chmixture[, "y.Mean.6"], auto.layout=FALSE, ask=FALSE, 
                col="darkblue", main="E(Tooth 16)") 
  dev.off()
}


###################################################
### chunk number 40: Model0: Autocorrelation plots
###################################################
#line 889 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob15.ps", sep=""), width=7, height=10, 
             horizontal=FALSE)
  par(mfrow=c(2, 3), bty="n")
  autocorr.plot(chmixture[, "y.SD.1"], auto.layout=FALSE, ask=FALSE, 
                col="darkgreen", main="SD(Tooth 11)")
  autocorr.plot(chmixture[, "y.SD.2"], auto.layout=FALSE, ask=FALSE, 
                col="darkgreen", main="SD(Tooth 12)")
  autocorr.plot(chmixture[, "y.SD.3"], auto.layout=FALSE, ask=FALSE, 
                col="darkgreen", main="SD(Tooth 13)")
  autocorr.plot(chmixture[, "y.SD.4"], auto.layout=FALSE, ask=FALSE, 
                col="darkgreen", main="SD(Tooth 14)")
  autocorr.plot(chmixture[, "y.SD.5"], auto.layout=FALSE, ask=FALSE, 
                col="darkgreen", main="SD(Tooth 15)")
  autocorr.plot(chmixture[, "y.SD.6"], auto.layout=FALSE, ask=FALSE, 
                col="darkgreen", main="SD(Tooth 16)")
  dev.off()
}


###################################################
### chunk number 41: Model0: Autocorrelation plots
###################################################
#line 927 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob16.ps", sep=""), width=10, height=7, 
             horizontal=FALSE)
  par(bty="n", mar=c(1, 1, 4, 0)+0.1)
  layout(matrix(c(1,2,3,4,5, 0,6,7,8,9, 0,0,10,11,12, 
                  0,0,0,13,14, 0,0,0,0,15), nrow=5, byrow=TRUE))
  autocorr.plot(chmixture[, "y.Corr.2.1"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 2-1")
  autocorr.plot(chmixture[, "y.Corr.3.1"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 3-1")
  autocorr.plot(chmixture[, "y.Corr.4.1"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 4-1")
  autocorr.plot(chmixture[, "y.Corr.5.1"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 5-1")
  autocorr.plot(chmixture[, "y.Corr.6.1"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 6-1")
  autocorr.plot(chmixture[, "y.Corr.3.2"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 3-2")
  autocorr.plot(chmixture[, "y.Corr.4.2"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 4-2")
  autocorr.plot(chmixture[, "y.Corr.5.2"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 5-2")
  autocorr.plot(chmixture[, "y.Corr.6.2"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 6-2")
  autocorr.plot(chmixture[, "y.Corr.4.3"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 4-3")
  autocorr.plot(chmixture[, "y.Corr.5.3"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 5-3")
  autocorr.plot(chmixture[, "y.Corr.6.3"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 6-3")
  autocorr.plot(chmixture[, "y.Corr.5.4"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 5-4")
  autocorr.plot(chmixture[, "y.Corr.6.4"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 6-4")
  autocorr.plot(chmixture[, "y.Corr.6.5"], auto.layout=FALSE, 
                ask=FALSE, col="darkgreen", main="Corr 6-5")
  dev.off()
}


###################################################
### chunk number 42: Model0: Autocorrelation plots
###################################################
#line 979 "Tandmob.Rnw"
if (RUN.ALLOUT){
  postscript(paste(FIGKEEPDIR, "figTandmob17.ps", sep=""), width=10, height=7, 
             horizontal=FALSE)
  par(bty="n")
  layout(matrix(c(1,1,2,2,3,3, 0,4,4,5,5,0), nrow=2, byrow=TRUE))
  autocorr.plot(chgammaInv[, "gammaInv1"], auto.layout=FALSE, 
                ask=FALSE, col="brown", main="gamma^{-1}")
  autocorr.plot(chdeviance[, "LogL0"], auto.layout=FALSE, 
                ask=FALSE, col="red", main="Log(L0)")
  autocorr.plot(chdeviance[, "LogL1"], auto.layout=FALSE, 
                ask=FALSE, col="red", main="Log(L1)")
  autocorr.plot(chdeviance[, "dev.complete"], auto.layout=FALSE, 
                ask=FALSE, col="red", main="D(complete)")
  autocorr.plot(chdeviance[, "dev.observed"], auto.layout=FALSE, 
                ask=FALSE, col="red", main="D(observed)")
  dev.off()
}


###################################################
### chunk number 43: Fixed K MCMC
###################################################
#line 1020 "Tandmob.Rnw"
Seed <- c(770328, 770328, 770328, 770328, 770328, 
          770328, 780830, 761014, 770328, 770328)

Keep <- c("iter", "nMCMC", "dim", "prior", "init", "RJMCMC", 
          "scale", "freqK", "propK", "DIC", "moves", 
          "summ.y.Mean", "summ.y.SDCorr", 
          "summ.z.Mean", "summ.z.SDCorr")

if (RUN.TIMECONSUMING.CODE){
  ModelK <- list()
  PDensUniModelK <- list()
  PDensBiModelK <- list()
  for (K in 1:10){  
    set.seed(Seed[K])      
    cat(paste("K = ", K, "\n-------------------------------\n", sep=""))
    cat("Seed is", Seed[K], "\n")  
    Prior0 <- list(priorK="fixed", Kmax=K)
    Model0 <- NMixMCMC(y0=y0, y1=y1, censor=censor, prior=Prior0, 
                       nMCMC=nMCMC, scale=list(shift=0, scale=1), PED=TRUE)
  
    cat(paste("\nComputation of marginal pred. densities started on ", date(), 
              "\n", sep=""))  
    PDensUniModelK[[k]] <- list()
    PDensUniModelK[[k]][[1]] <- NMixPredDensMarg(Model0[[1]], grid=ygridUni)
    PDensUniModelK[[k]][[2]] <- NMixPredDensMarg(Model0[[2]], grid=ygridUni)
    cat(paste("Computation of marginal pred. densities finished on ", date(), 
              "\n\n", sep=""))  

    cat(paste("Computation of joint pred. densities started on ", date(), 
              "\n", sep=""))    
    PDensBiModelK[[k]] <- list()
    PDensBiModelK[[k]][[1]] <- NMixPredDensJoint2(Model0[[1]], grid=ygridBi)
    PDensBiModelK[[k]][[2]] <- NMixPredDensJoint2(Model0[[2]], grid=ygridBi)
    cat(paste("Computation of joint pred. densities finished on ", date(), 
              "\n\n\n", sep=""))      

    ModelK[[k]][[1]] <- Model0[[1]][Keep]
    ModelK[[k]][[2]] <- Model0[[2]][Keep]  
    class(ModelK[[k]][[1]]) <- class(ModelK[[k]][[2]]) <- "NMixMCMC"

    rm(list="Model0")
  }  
}


###################################################
### chunk number 44: Fixed K: PED and DIC
###################################################
#line 1068 "Tandmob.Rnw"
PED <- ModelK[[1]]$PED
DIC <- list(Chain1=ModelK[[1]][[1]]$DIC, Chain2=ModelK[[1]][[2]]$DIC)
for (k in 2:length(ModelK)){
  PED <- rbind(PED, ModelK[[k]]$PED)
  DIC[[1]] <- rbind(DIC[[1]], ModelK[[k]][[1]]$DIC)
  DIC[[2]] <- rbind(DIC[[2]], ModelK[[k]][[2]]$DIC)  
}
rownames(PED) <- paste("K=", 1:length(ModelK), sep="")  
rownames(DIC[[1]]) <- rownames(DIC[[2]]) <- paste("K=", 1:length(ModelK), sep="")  


###################################################
### chunk number 45: Fixed K: print PED
###################################################
#line 1079 "Tandmob.Rnw"
print(PED)


###################################################
### chunk number 46: Fixed K: print DIC
###################################################
#line 1082 "Tandmob.Rnw"
print(DIC)


###################################################
### chunk number 47: FixedK: Plots of univariate predictive densities for different values of K
###################################################
#line 1089 "Tandmob.Rnw"
if (RUN.ALLOUT){
  CH <- 1
  for (fig in 1:3){
    postscript(paste(FIGKEEPDIR, "figTandmob", 17+fig, ".ps", sep=""), 
               width=7, height=10, horizontal=FALSE)  
    par(mfrow=c(2, 1), bty="n")
    for (tt in (2*(fig-1)+1):(2*fig)){
      plot(PDensUniModelK[[Kshow]][[CH]]$x[[tt]], 
           PDensUniModelK[[Kshow]][[CH]]$dens[[tt]], 
           type="n", xlab="Age (years)", ylab="Density", 
           main=paste("Tooth ", 10+tt, sep=""))    
      for (k in 1:10){
        lines(PDensUniModelK[[k]][[CH]]$x[[tt]], 
              PDensUniModelK[[k]][[CH]]$dens[[tt]], 
              col="blue") 
      }  
      lines(PDensUniModelK[[Kshow]][[CH]]$x[[tt]], 
            PDensUniModelK[[Kshow]][[CH]]$dens[[tt]], 
            col="red", lwd=2)     
    }
    dev.off()
  }  
}


###################################################
### chunk number 48: Save results
###################################################
#line 1137 "Tandmob.Rnw"
if (RUN.TIMECONSUMING.CODE){
  save(list="Model0", 
       file=paste(RESULTDIR, "/Tandmob-Model0", Kshow, ".RData", sep=""))            

  PDensBiModelK <- list()
  PDensBiModelK[[Kshow]] <- PDensBiModel0
  save(list="PDensBiModelK", 
       file=paste(RESULT2DIR, "/Tandmob-PDensBiModelK0", Kshow, ".RData", sep=""))  
  
  save(list=c("ModelK", "PDensUniModelK", "PDensBiModelK"), 
       file=paste(RESULT2DIR, "/Tandmob-Result.RData", sep=""))      
}  


