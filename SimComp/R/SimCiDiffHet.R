SimCiDiffHet <-
function(trlist, grp, ntr, nep, ssvec, Cmat, alternative, conf.level) {


ncomp <- nrow(Cmat)                                                 # number of comparisons

meanmat <- varmat <- matrix(nrow=ntr, ncol=nep)
for (i in 1:ntr) { for (j in 1:nep) {
  meanmat[i,j]=mean(trlist[[i]][,j]); varmat[i,j]=var(trlist[[i]][,j]) }}
estimate <- Cmat%*%meanmat

defrmat <- matrix(nrow=ncomp, ncol=nep)
for (j in 1:nep) { for (z in 1:ncomp) {
defrmat[z,j]=( (sum((Cmat[z,])^2*varmat[,j]/ssvec))^2 ) / 
             sum( ( (Cmat[z,])^4*varmat[,j]^2 ) / ( ssvec^2*(ssvec-1) ) ) }}
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
      Rpart[i,h]=( t(Cmat[z,])%*%diag(unlist( lapply( X=CovMatDat,FUN=function(x){x[i,h]} ) ))%*%M%*%(Cmat[w,]) ) /
                 sqrt( ( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%(Cmat[z,]) ) *
                       ( t(Cmat[w,])%*%diag(varmat[,h])%*%M%*%(Cmat[w,]) ) ) }
    }
    Rrow <- cbind(Rrow,Rpart)
  }
  R <- rbind(R, Rrow)                                               # correlation matrix for test.stat
}
diag(R) <- 1

lower <- upper <- lower.raw <- upper.raw <- matrix(nrow=ncomp,ncol=nep)

if (alternative=="greater") {
  for (z in 1:ncomp) { for (i in 1:nep) {
    lo1malqu <- qmvt(conf.level,tail="lower.tail",df=as.integer(defrvec[z]),corr=R)$quantile
    univarqu <- qt(p=conf.level, df=defrmat[z,i])
    upper[z,i] <- upper.raw[z,i] <- Inf
    lower[z,i]     <- t(Cmat[z,])%*%meanmat[,i] - lo1malqu * sqrt( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%Cmat[z,] )
    lower.raw[z,i] <- t(Cmat[z,])%*%meanmat[,i] - univarqu * sqrt( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%Cmat[z,] )
  }}
}
if (alternative=="less") {
  for (z in 1:ncomp) { for (i in 1:nep) {
    up1malqu <- qmvt(conf.level,tail="upper.tail",df=as.integer(defrvec[z]),corr=R)$quantile
    univarqu <- qt(p=1-conf.level, df=defrmat[z,i])
    upper[z,i]     <- t(Cmat[z,])%*%meanmat[,i] - up1malqu * sqrt( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%Cmat[z,] )
    upper.raw[z,i] <- t(Cmat[z,])%*%meanmat[,i] - univarqu * sqrt( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%Cmat[z,] )
    lower[z,i] <- lower.raw[z,i] <- -Inf
  }}
}
if (alternative=="two.sided") {
  for (z in 1:ncomp) { for (i in 1:nep) {
    ts1malqu <- qmvt(conf.level,tail="both.tails",df=as.integer(defrvec[z]),corr=R)$quantile
    univarqu <- qt(p=1-(1-conf.level)/2, df=defrmat[z,i])
    upper[z,i]     <- t(Cmat[z,])%*%meanmat[,i] + ts1malqu * sqrt( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%Cmat[z,] )
    upper.raw[z,i] <- t(Cmat[z,])%*%meanmat[,i] + univarqu * sqrt( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%Cmat[z,] )
    lower[z,i]     <- t(Cmat[z,])%*%meanmat[,i] - ts1malqu * sqrt( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%Cmat[z,] )
    lower.raw[z,i] <- t(Cmat[z,])%*%meanmat[,i] - univarqu * sqrt( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%Cmat[z,] )
  }}
}

list(estimate=estimate, lower.raw=lower.raw, upper.raw=upper.raw, lower=lower, upper=upper,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defrvec, 
     Cmat=Cmat, alternative=alternative, conf.level=conf.level)


}
