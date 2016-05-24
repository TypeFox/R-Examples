SimCiRatHom <-
function(trlist, grp, ntr, nep, ssvec, Num.Contrast, Den.Contrast, alternative, conf.level) {


CMat <- Num.Contrast
DMat <- Den.Contrast
ncomp <- nrow(CMat)                                                 # number of comparisons

meanmat <- matrix(nrow=ntr, ncol=nep)
for (i in 1:ntr) { for (j in 1:nep) {
  meanmat[i,j]=mean(trlist[[i]][,j]) }}
if (any(meanmat<0)) {
  cat("Warning: At least one sample mean is negative; check whether the test direction", "\n",
      "is still correct", "\n")
}
estimate <- CMat%*%meanmat/(DMat%*%meanmat)

defr <- sum(ssvec)-ntr                                              # degrees of freedom

CovMatDat <- matrix(rep(0,nep*nep),nrow=nep)                        # common covariance matrix of the data
for (i in 1:ntr) { CovMatDat <- CovMatDat+(ssvec[i]-1)*cov(trlist[[i]]) }
CovMatDat <- CovMatDat/defr
CorrMatDat <- cov2cor(CovMatDat)                                    # common correlation matrix of the data

M <- diag(1/ssvec)
R <- NULL
for (z in 1:ncomp) { 
  Rrow <- NULL
  for (w in 1:ncomp) {
    Rpart <- matrix(nrow=nep,ncol=nep)
    for (i in 1:nep) { for (h in 1:nep) {
      Rpart[i,h]=CorrMatDat[i,h] * ( t(CMat[z,]-estimate[z,i]*DMat[z,])%*%M%*%(CMat[w,]-estimate[w,h]*DMat[w,]) ) /
                 sqrt( (t(CMat[z,]-estimate[z,i]*DMat[z,])%*%M%*%(CMat[z,]-estimate[z,i]*DMat[z,])) * 
                       (t(CMat[w,]-estimate[w,h]*DMat[w,])%*%M%*%(CMat[w,]-estimate[w,h]*DMat[w,])) ) }
    }
    Rrow <- cbind(Rrow,Rpart)
  }
  R <- rbind(R, Rrow)                                               # correlation matrix for test.stat
}
diag(R) <- 1

Azi     <- Bzi     <- Czi     <- Discrimi     <- lower     <- upper     <- matrix(nrow=ncomp,ncol=nep)
Azi.raw <- Bzi.raw <- Czi.raw <- Discrimi.raw <- lower.raw <- upper.raw <- matrix(nrow=ncomp,ncol=nep)
NSD <- 0

if (alternative=="greater") {
  lo1malqu <- qmvt(conf.level,tail="lower.tail",df=defr,corr=R)$quantile
  univarqu <- qt(p=conf.level, df=defr)
  for (z in 1:ncomp) { for (i in 1:nep) {
    Azi[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - lo1malqu^2 * diag(CovMatDat)[i] * ( t(DMat[z,])%*%M%*%DMat[z,] )
    Bzi[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                - lo1malqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%DMat[z,] ) )
    Czi[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - lo1malqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%CMat[z,] )
    Discrimi[z,i] <- Bzi[z,i]^2 - 4*Azi[z,i]*Czi[z,i]
    if ( (Azi[z,i]>0) & (Discrimi[z,i]>=0) ) {
      upper[z,i] <- Inf
      lower[z,i] <- (-Bzi[z,i]-sqrt(Discrimi[z,i])) / (2*Azi[z,i])
    } else {
      upper[z,i] <- Inf; lower[z,i] <- -Inf; NSD <- NSD+1
    }
    Azi.raw[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * diag(CovMatDat)[i] * ( t(DMat[z,])%*%M%*%DMat[z,] )
    Bzi.raw[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                    - univarqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%DMat[z,] ) )
    Czi.raw[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%CMat[z,] )
    Discrimi.raw[z,i] <- Bzi.raw[z,i]^2 - 4*Azi.raw[z,i]*Czi.raw[z,i]
    if ( (Azi.raw[z,i]>0) & (Discrimi.raw[z,i]>=0) ) {
      upper.raw[z,i] <- Inf
      lower.raw[z,i] <- (-Bzi.raw[z,i]-sqrt(Discrimi.raw[z,i])) / (2*Azi.raw[z,i])
    } else {
      upper.raw[z,i] <- Inf; lower.raw[z,i] <- -Inf
    }
  }}
}
if (alternative=="less") {
  up1malqu <- qmvt(conf.level,tail="upper.tail",df=defr,corr=R)$quantile
  univarqu <- qt(p=1-conf.level, df=defr)
  for (z in 1:ncomp) { for (i in 1:nep) {
    Azi[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - up1malqu^2 * diag(CovMatDat)[i] * ( t(DMat[z,])%*%M%*%DMat[z,] )
    Bzi[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                - up1malqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%DMat[z,] ) )
    Czi[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - up1malqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%CMat[z,] )
    Discrimi[z,i] <- Bzi[z,i]^2 - 4*Azi[z,i]*Czi[z,i]
    if ( (Azi[z,i]>0) & (Discrimi[z,i]>=0) ) {
      upper[z,i] <- (-Bzi[z,i]+sqrt(Discrimi[z,i])) / (2*Azi[z,i])
      lower[z,i] <- -Inf
    } else {
      upper[z,i] <- Inf; lower[z,i] <- -Inf; NSD <- NSD+1
    }
    Azi.raw[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * diag(CovMatDat)[i] * ( t(DMat[z,])%*%M%*%DMat[z,] )
    Bzi.raw[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                    - univarqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%DMat[z,] ) )
    Czi.raw[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%CMat[z,] )
    Discrimi.raw[z,i] <- Bzi.raw[z,i]^2 - 4*Azi.raw[z,i]*Czi.raw[z,i]
    if ( (Azi.raw[z,i]>0) & (Discrimi.raw[z,i]>=0) ) {
      upper.raw[z,i] <- (-Bzi.raw[z,i]+sqrt(Discrimi.raw[z,i])) / (2*Azi.raw[z,i])
      lower.raw[z,i] <- -Inf
    } else {
      upper.raw[z,i] <- Inf; lower.raw[z,i] <- -Inf
    }
  }}
}
if (alternative=="two.sided") {
  ts1malqu <- qmvt(conf.level,tail="both.tails",df=defr,corr=R)$quantile
  univarqu <- qt(p=1-(1-conf.level)/2, df=defr)
  for (z in 1:ncomp) { for (i in 1:nep) {
    Azi[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - ts1malqu^2 * diag(CovMatDat)[i] * ( t(DMat[z,])%*%M%*%DMat[z,] )
    Bzi[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                - ts1malqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%DMat[z,] ) )
    Czi[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - ts1malqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%CMat[z,] )
    Discrimi[z,i] <- Bzi[z,i]^2 - 4*Azi[z,i]*Czi[z,i]
    if ( (Azi[z,i]>0) & (Discrimi[z,i]>=0) ) {
      upper[z,i] <- (-Bzi[z,i]+sqrt(Discrimi[z,i])) / (2*Azi[z,i])
      lower[z,i] <- (-Bzi[z,i]-sqrt(Discrimi[z,i])) / (2*Azi[z,i])
    } else {
      upper[z,i] <- Inf; lower[z,i] <- -Inf; NSD <- NSD+1
    }
    Azi.raw[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * diag(CovMatDat)[i] * ( t(DMat[z,])%*%M%*%DMat[z,] )
    Bzi.raw[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                    - univarqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%DMat[z,] ) )
    Czi.raw[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * diag(CovMatDat)[i] * ( t(CMat[z,])%*%M%*%CMat[z,] )
    Discrimi.raw[z,i] <- Bzi.raw[z,i]^2 - 4*Azi.raw[z,i]*Czi.raw[z,i]
    if ( (Azi.raw[z,i]>0) & (Discrimi.raw[z,i]>=0) ) {
      upper.raw[z,i] <- (-Bzi.raw[z,i]+sqrt(Discrimi.raw[z,i])) / (2*Azi.raw[z,i])
      lower.raw[z,i] <- (-Bzi.raw[z,i]-sqrt(Discrimi.raw[z,i])) / (2*Azi.raw[z,i])
    } else {
      upper.raw[z,i] <- Inf; lower.raw[z,i] <- -Inf
    }
  }}
}

list(estimate=estimate, NSD=NSD, lower.raw=lower.raw, upper.raw=upper.raw, lower=lower, upper=upper,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr, 
     Num.Contrast=CMat, Den.Contrast=DMat, alternative=alternative, conf.level=conf.level)


}
