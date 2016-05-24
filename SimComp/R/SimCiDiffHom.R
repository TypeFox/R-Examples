SimCiDiffHom <-
function(trlist, grp, ntr, nep, ssvec, Cmat, alternative, conf.level) {


ncomp <- nrow(Cmat)                                                 # number of comparisons

meanmat <- matrix(nrow=ntr, ncol=nep)
for (i in 1:ntr) { for (j in 1:nep) {
  meanmat[i,j]=mean(trlist[[i]][,j]) }}
estimate <- Cmat%*%meanmat

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
      Rpart[i,h]=CorrMatDat[i,h] * (t(Cmat[z,])%*%M%*%Cmat[w,]) /
                 sqrt( (t(Cmat[z,])%*%M%*%Cmat[z,]) * ((Cmat[w,])%*%M%*%Cmat[w,]) ) }
    }
    Rrow <- cbind(Rrow,Rpart)
  }
  R <- rbind(R, Rrow)                                               # correlation matrix for test.stat
}
diag(R) <- 1

lower <- upper <- lower.raw <- upper.raw <- matrix(nrow=ncomp,ncol=nep)

if (alternative=="greater") {
  lo1malqu <- qmvt(conf.level,tail="lower.tail",df=defr,corr=R)$quantile
  univarqu <- qt(p=conf.level, df=defr)
  for (z in 1:ncomp) { for (i in 1:nep) {
    upper[z,i] <- upper.raw[z,i] <- Inf
    lower[z,i]     <- t(Cmat[z,])%*%meanmat[,i] - lo1malqu * sqrt( diag(CovMatDat)[i]*( t(Cmat[z,])%*%M%*%Cmat[z,] ) )
    lower.raw[z,i] <- t(Cmat[z,])%*%meanmat[,i] - univarqu * sqrt( diag(CovMatDat)[i]*( t(Cmat[z,])%*%M%*%Cmat[z,] ) )
  }}
}
if (alternative=="less") {
  up1malqu <- qmvt(conf.level,tail="upper.tail",df=defr,corr=R)$quantile
  univarqu <- qt(p=1-conf.level, df=defr)
  for (z in 1:ncomp) { for (i in 1:nep) {
    upper[z,i]     <- t(Cmat[z,])%*%meanmat[,i] - up1malqu * sqrt( diag(CovMatDat)[i]*( t(Cmat[z,])%*%M%*%Cmat[z,] ) )
    upper.raw[z,i] <- t(Cmat[z,])%*%meanmat[,i] - univarqu * sqrt( diag(CovMatDat)[i]*( t(Cmat[z,])%*%M%*%Cmat[z,] ) )
    lower[z,i] <- lower.raw[z,i] <- -Inf
  }}
}
if (alternative=="two.sided") {
  ts1malqu <- qmvt(conf.level,tail="both.tails",df=defr,corr=R)$quantile
  univarqu <- qt(p=1-(1-conf.level)/2, df=defr)
  for (z in 1:ncomp) { for (i in 1:nep) {
    upper[z,i]     <- t(Cmat[z,])%*%meanmat[,i] + ts1malqu * sqrt( diag(CovMatDat)[i]*( t(Cmat[z,])%*%M%*%Cmat[z,] ) )
    upper.raw[z,i] <- t(Cmat[z,])%*%meanmat[,i] + univarqu * sqrt( diag(CovMatDat)[i]*( t(Cmat[z,])%*%M%*%Cmat[z,] ) )
    lower[z,i]     <- t(Cmat[z,])%*%meanmat[,i] - ts1malqu * sqrt( diag(CovMatDat)[i]*( t(Cmat[z,])%*%M%*%Cmat[z,] ) )
    lower.raw[z,i] <- t(Cmat[z,])%*%meanmat[,i] - univarqu * sqrt( diag(CovMatDat)[i]*( t(Cmat[z,])%*%M%*%Cmat[z,] ) )
  }}
}

list(estimate=estimate, lower.raw=lower.raw, upper.raw=upper.raw, lower=lower, upper=upper,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr, 
     Cmat=Cmat, alternative=alternative, conf.level=conf.level)


}
