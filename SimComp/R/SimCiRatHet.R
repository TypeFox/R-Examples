SimCiRatHet <-
function(trlist, grp, ntr, nep, ssvec, Num.Contrast, Den.Contrast, alternative, conf.level) {


CMat <- Num.Contrast
DMat <- Den.Contrast
ncomp <- nrow(CMat)                                                 # number of comparisons

meanmat <- varmat <- matrix(nrow=ntr, ncol=nep)
for (i in 1:ntr) { for (j in 1:nep) {
  meanmat[i,j]=mean(trlist[[i]][,j]); varmat[i,j]=var(trlist[[i]][,j]) }}
if (any(meanmat<0)) {
  cat("Warning: At least one sample mean is negative; check whether the test direction", "\n",
      "is still correct", "\n")
}
estimate <- CMat%*%meanmat/(DMat%*%meanmat)

defrmat <- matrix(nrow=ncomp, ncol=nep)
for (j in 1:nep) { for (z in 1:ncomp) {
defrmat[z,j]=( (sum((CMat[z,]-estimate[z,j]*DMat[z,])^2*varmat[,j]/ssvec))^2 ) / 
             sum( ( (CMat[z,]-estimate[z,j]*DMat[z,])^4*varmat[,j]^2 ) / ( ssvec^2*(ssvec-1) ) ) }}
defrmat[defrmat<2] <- 2                                             # to be well-defined
defrvec <- apply(X=defrmat, MARGIN=1, FUN=min)                      # minimum over the rows/endpoints

CovMatDat <- CorrMatDat <- list()                                   # list of covariance/correlation matrices of the data
for (i in 1:ntr) { CovMatDat[[i]]  <- cov(trlist[[i]])
                   CorrMatDat[[i]] <- cov2cor(CovMatDat[[i]]) }

M <- diag(1/ssvec)
R <- NULL
for (z in 1:ncomp) {
  Rrow <- NULL
  for (w in 1:ncomp) {
    Rpart <- matrix(nrow=nep,ncol=nep)
    for (i in 1:nep) { for (h in 1:nep) {
      Rpart[i,h]=( t(CMat[z,]-estimate[z,i]*DMat[z,])%*%
                   diag(unlist( lapply( X=CovMatDat,FUN=function(x){x[i,h]} ) ))%*%M%*%(CMat[w,]-estimate[w,h]*DMat[w,]) ) /
                 sqrt( ( t(CMat[z,]-estimate[z,i]*DMat[z,])%*%diag(varmat[,i])%*%M%*%(CMat[z,]-estimate[z,i]*DMat[z,]) ) *
                       ( t(CMat[w,]-estimate[w,h]*DMat[w,])%*%diag(varmat[,h])%*%M%*%(CMat[w,]-estimate[w,h]*DMat[w,]) ) ) }
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
  for (z in 1:ncomp) { for (i in 1:nep) {
    lo1malqu <- qmvt(conf.level,tail="lower.tail",df=as.integer(defrvec[z]),corr=R)$quantile
    univarqu <- qt(p=conf.level, df=defrmat[z,i])
    Azi[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - lo1malqu^2 * ( t(DMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] )
    Bzi[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                - lo1malqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] ) )
    Czi[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - lo1malqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%CMat[z,] )
    Discrimi[z,i] <- Bzi[z,i]^2 - 4*Azi[z,i]*Czi[z,i]
    if ( (Azi[z,i]>0) & (Discrimi[z,i]>=0) ) {
      upper[z,i] <- Inf
      lower[z,i] <- (-Bzi[z,i]-sqrt(Discrimi[z,i])) / (2*Azi[z,i])
    } else {
      upper[z,i] <- Inf; lower[z,i] <- -Inf; NSD <- NSD+1
    }
    Azi.raw[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * ( t(DMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] )
    Bzi.raw[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                    - univarqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] ) )
    Czi.raw[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%CMat[z,] )
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
  for (z in 1:ncomp) { for (i in 1:nep) {
    up1malqu <- qmvt(conf.level,tail="upper.tail",df=as.integer(defrvec[z]),corr=R)$quantile
    univarqu <- qt(p=1-conf.level, df=defrmat[z,i])
    Azi[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - up1malqu^2 * ( t(DMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] )
    Bzi[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                - up1malqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] ) )
    Czi[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - up1malqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%CMat[z,] )
    Discrimi[z,i] <- Bzi[z,i]^2 - 4*Azi[z,i]*Czi[z,i]
    if ( (Azi[z,i]>0) & (Discrimi[z,i]>=0) ) {
      upper[z,i] <- (-Bzi[z,i]+sqrt(Discrimi[z,i])) / (2*Azi[z,i])
      lower[z,i] <- -Inf
    } else {
      upper[z,i] <- Inf; lower[z,i] <- -Inf; NSD <- NSD+1
    }
    Azi.raw[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * ( t(DMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] )
    Bzi.raw[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                    - univarqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] ) )
    Czi.raw[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%CMat[z,] )
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
  for (z in 1:ncomp) { for (i in 1:nep) {
    ts1malqu <- qmvt(conf.level,tail="both.tails",df=as.integer(defrvec[z]),corr=R)$quantile
    univarqu <- qt(p=1-(1-conf.level)/2, df=defrmat[z,i])
    Azi[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - ts1malqu^2 * ( t(DMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] )
    Bzi[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                - ts1malqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] ) )
    Czi[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - ts1malqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%CMat[z,] )
    Discrimi[z,i] <- Bzi[z,i]^2 - 4*Azi[z,i]*Czi[z,i]
    if ( (Azi[z,i]>0) & (Discrimi[z,i]>=0) ) {
      upper[z,i] <- (-Bzi[z,i]+sqrt(Discrimi[z,i])) / (2*Azi[z,i])
      lower[z,i] <- (-Bzi[z,i]-sqrt(Discrimi[z,i])) / (2*Azi[z,i])
    } else {
      upper[z,i] <- Inf; lower[z,i] <- -Inf; NSD <- NSD+1
    }
    Azi.raw[z,i] <- ( t(DMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * ( t(DMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] )
    Bzi.raw[z,i] <- -2 * ( (t(CMat[z,])%*%meanmat[,i]) * (t(DMat[z,])%*%meanmat[,i])
                                                    - univarqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%DMat[z,] ) )
    Czi.raw[z,i] <- ( t(CMat[z,])%*%meanmat[,i] )^2 - univarqu^2 * ( t(CMat[z,])%*%diag(varmat[,i])%*%M%*%CMat[z,] )
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
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defrvec, 
     Num.Contrast=CMat, Den.Contrast=DMat, alternative=alternative, conf.level=conf.level)


}
