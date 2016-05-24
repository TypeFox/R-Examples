SimTestDiffHet <-
function(trlist, grp, ntr, nep, ssvec, Cmat, alternative, Margin) {


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

test.stat <- p.val.adj <- p.val.raw <- matrix(nrow=ncomp, ncol=nep) # matrices of test statistics and p.vals
for (z in 1:ncomp) { for (i in 1:nep) {
  test.stat[z,i]=( (t(Cmat[z,])%*%meanmat[,i]) - Margin[z,i] ) /
                 sqrt( t(Cmat[z,])%*%diag(varmat[,i])%*%M%*%(Cmat[z,]) )
  if (alternative=="greater") {
    p.val.adj[z,i]=1-pmvt(lower=-Inf,upper=rep(test.stat[z,i],times=ncomp*nep),df=as.integer(defrvec[z]),corr=R)[1]
    p.val.raw[z,i]=pt(q=test.stat[z,i],df=defrmat[z,i],lower.tail=FALSE) }
  if (alternative=="less") {
    p.val.adj[z,i]=1-pmvt(lower=rep(test.stat[z,i],times=ncomp*nep),upper=Inf,df=as.integer(defrvec[z]),corr=R)[1]
    p.val.raw[z,i]=pt(q=test.stat[z,i],df=defrmat[z,i],lower.tail=TRUE) }
  if (alternative=="two.sided") {
    p.val.adj[z,i]=1-pmvt(lower=rep(-abs(test.stat[z,i]),times=ncomp*nep),
                   upper=rep(abs(test.stat[z,i]),times=ncomp*nep),df=as.integer(defrvec[z]),corr=R)[1]
    p.val.raw[z,i]=min(pt(q=abs(test.stat[z,i]),df=defrmat[z,i],lower.tail=FALSE)*2,1) }
}}

list(estimate=estimate, statistic=test.stat, p.val.raw=p.val.raw, p.val.adj=p.val.adj,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defrvec,
     Cmat=Cmat, Margin=Margin, alternative=alternative)


}
