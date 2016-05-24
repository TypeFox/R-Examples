ACE <- function(mzDat=mzData,dzDat=dzData,type="raw",selV=selVars)
{
  require(OpenMx)
  twinACEModel <- mxModel("ACE",
  mxMatrix("Full", 1, 1, TRUE, .6, "a", name="X"),
  mxMatrix("Full", 1, 1, TRUE, .6, "c", name="Y"),
  mxMatrix("Full", 1, 1, TRUE, .6, "e", name="Z"),
  mxAlgebra(X %*% t(X), "A"),
  mxAlgebra(Y %*% t(Y), "C"),
  mxAlgebra(Z %*% t(Z), "E"),
  mxAlgebra(A+C+E, name="V"),
  mxMatrix("Full", 1, 2, TRUE, 20, "mean", name="expMean"),
  mxAlgebra(rbind(cbind(A+C+E, A+C), cbind(A+C, A+C+E)), "expCovMZ"),
  mxAlgebra(rbind(cbind(A+C+E, 0.5%x%A+C), cbind(0.5%x%A+C, A+C+E)), "expCovDZ"),
  mxModel("MZ", mxData(mzDat, type), mxFIMLObjective("ACE.expCovMZ", "ACE.expMean", selV)),
  mxModel("DZ", mxData(dzDat, type), mxFIMLObjective("ACE.expCovDZ", "ACE.expMean", selV)),
  mxAlgebra(MZ.objective + DZ.objective, name="twin"),
  mxAlgebraObjective("twin"))

  twinACEFit <- mxRun(twinACEModel, silent=TRUE)
  exp_ACE <- mxEval(rbind(expCovMZ,expCovDZ,expMean), twinACEFit)
  est_ACE <- mxEval(cbind(A,C,E,A/V,C/V,E/V), twinACEFit)
  LL_ACE <- mxEval(objective, twinACEFit)
  rownames(exp_ACE) <- c('CovMZT1','CovMZT2','CovDZT1','CovDZT2','Mean')
  colnames(exp_ACE) <- c('T1','T2')
  rownames(est_ACE) <- 'ACE'
  colnames(est_ACE) <- c('a','c','e','a^2','c^2','e^2')

  twinAEModel <- mxModel(twinACEModel, mxMatrix("Full", 1, 1, FALSE, 0, "c", name="Y"))
  twinAEFit <- mxRun(twinAEModel, silent=TRUE)
  exp_AE <- mxEval(rbind(expCovMZ,expCovDZ,expMean), twinAEFit)
  est_AE <- mxEval(cbind(A,C,E,A/V,C/V,E/V), twinAEFit)
  rownames(est_AE) <- 'AE'
  LL_AE <- mxEval(objective, twinAEFit)
  LRT_ACE_AE <- LL_AE - LL_ACE

  twinCEModel <- mxModel(twinACEModel, mxMatrix("Full", 1, 1, FALSE, 0, "a", name="X"))
  twinCEFit <- mxRun(twinCEModel, silent=TRUE)
  exp_CE <- mxEval(rbind(expCovMZ,expCovDZ,expMean), twinCEFit)
  est_CE <- mxEval(cbind(A,C,E,A/V,C/V,E/V), twinCEFit)
  rownames(est_CE) <- 'CE'
  LL_CE <- mxEval(objective, twinCEFit)
  LRT_ACE_CE <- LL_CE - LL_ACE

  twinEModel <- mxModel(twinAEModel, mxMatrix("Full", 1, 1, FALSE, 0, "a", name="X"))
  twinEFit <- mxRun(twinEModel, silent=TRUE)
  exp_E <- mxEval(rbind(expCovMZ,expCovDZ,expMean), twinEFit)
  est_E <- mxEval(cbind(A,C,E,A/V,C/V,E/V), twinEFit)
  rownames(est_E) <- 'E'
  LL_E <- mxEval(objective, twinEFit)
  LRT_ACE_E <- LL_E - LL_ACE

  exp <- cbind(exp_ACE,exp_AE,exp_CE,exp_E)
  est <- rbind(est_ACE,est_AE,est_CE,est_E)
  lls <- rbind(cbind(LL_ACE,0),cbind(LL_AE,LRT_ACE_AE),cbind(LL_CE,LRT_ACE_CE),cbind(LL_E,LRT_ACE_E))
  df <- c(NA,1,1,2)
  lls <- cbind(lls,pchisq(2*lls[,2],df,lower.tail=FALSE))
  rownames(lls) <- c("l(ACE)","l(AE),lrt(AE,ACE)","l(CE),lrt(CE,ACE)","l(E),lrt(ACE,E)")
  invisible(list(exp=exp,est=est,lls=lls))
}

# Modified from code by Hermine Maes on 28/4/2010 MRC-Epid JHZ

k <- function(r,N,adjust=TRUE)
{
  r2 <- r^2
  n <- N-1
  k1 <- ifelse(adjust,r-r*(1-r2)/2/n,r)
  k2 <- (1-r2)^2/n*(1+11*r2/2/n)
  invisible(c(k1,k2))
}

# Implemented according to Keeping ES. Introduction to Statistical Inference, Dover Pulications, Inc. 1995

ACEr <- function(mzDat=mzData,dzDat=dzData,selV=selVars)
{
  rmz <- cor(mzDat[selV[1]],mzDat[selV[2]], use="complete")
  rdz <- cor(dzDat[selV[1]],dzDat[selV[2]], use="complete")
  nmz <- length(!is.na(c(mzDat[selV[1]],mzDat[selV[2]])))
  ndz <- length(!is.na(c(dzDat[selV[1]],dzDat[selV[2]])))
  kmz <- k(rmz,nmz)
  k1mz <- kmz[1]
  k2mz <- kmz[2]
  kdz <- k(rdz,ndz)
  k1dz <- kdz[1]
  k2dz <- kdz[2]
  h2 <- 2 * (k1mz - k1dz)
  vh <- 4 * (k2mz + k2dz)
  c2 <- 2 * k1dz - k1mz
  vc <- 4 * k2dz + k2mz
  e2 <- 1 - k1mz
  ve <- k2mz
  ACEr_est <- as.matrix(c(h2,c2,e2,vh,vc,ve))
  rownames(ACEr_est) <- c("h2","c2","e2","vh","vc","ve")
  invisible(t(ACEr_est))
}

# implemented according to wikipedia on 28/4/2010 MRC-Epid JHZ

ACE_CI <- function(mzData,dzData,n.sim=5,selV=selVars,verbose=TRUE)
{
options(echo=FALSE)

ACE_twinData <- ACE()
print(ACE_twinData)
ACEr_twinData <- ACEr(mzData,dzData,selV)
print(ACEr_twinData)

nmz <- dim(mzData)[1]
ndz <- dim(dzData)[1]
a <- ar <- vector()
set.seed(12345)
for(i in 1:n.sim)
{
  cat("\rRunning # ",i,"/", n.sim,"\r",sep="")
  sampled_mz <- sample(1:nmz, replace=TRUE)
  sampled_dz <- sample(1:ndz, replace=TRUE)
  mzDat <- mzData[sampled_mz,]
  dzDat <- dzData[sampled_dz,]
  ACE_i <- ACE(mzDat,dzDat,"raw",selV)
  if(verbose) print(ACE_i$est)
  a <- rbind(a,ACE_i$est[,4])
  ACEr_i <- ACEr(mzDat,dzDat,selV)
  if(verbose) print(ACEr_i)
  ar <- rbind(ar,ACEr_i)
}
a <- as.data.frame(a)
m <- mean(a, na.rm=TRUE)
s <- sd(a, na.rm=TRUE)
all <- data.frame(mean=m, sd=s, lcl=m-1.96*s, ucl=m+1.96*s)
cat("\nheritability and 95CI estimates by models\n\n")
print(all)
cat("\n\nheritability according to correlations\n\n")
ar <- as.data.frame(ar)
m <- mean(ar,na.rm=TRUE)
s <- sd(ar,na.rm=TRUE)
allr <- data.frame(mean=m,sd=s,lcl=m-1.96*s,ucl=m+1.96*s)
print(allr)
options(echo=TRUE)
}

# OpenMx documentation example
selVars <- c('bmi1','bmi2')

data(twinData, package="OpenMx")
mzData <- as.data.frame(subset(twinData, zyg==1, c(bmi1,bmi2)))
dzData <- as.data.frame(subset(twinData, zyg==3, c(bmi1,bmi2)))
ACE_CI(mzData,dzData)

# simulated data
library(mvtnorm)
n.sim <- 500
cat ("\nThe first study\n\n")
mzm <- as.data.frame(rmvnorm(195, c(22.75,22.75), matrix(2.66^2*c(1, 0.67, 0.67, 1), 2)))
dzm <- as.data.frame(rmvnorm(130, c(23.44,23.44), matrix(2.75^2*c(1, 0.32, 0.32, 1), 2)))
mzw <- as.data.frame(rmvnorm(384, c(21.44,21.44), matrix(3.08^2*c(1, 0.72, 0.72, 1), 2)))
dzw <- as.data.frame(rmvnorm(243, c(21.72,21.72), matrix(3.12^2*c(1, 0.33, 0.33, 1), 2)))
names(mzm) <- names(dzm) <- names(mzw) <- names(dzw) <- c("bmi1","bmi2")
ACE_CI(mzm,dzm,n.sim,verbose=FALSE)
ACE_CI(mzw,dzw,n.sim,verbose=FALSE)
