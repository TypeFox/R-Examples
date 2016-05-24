SimTestRatHom <-
function(trlist, grp, ntr, nep, ssvec, Num.Contrast, Den.Contrast, alternative, Margin) {


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
      Rpart[i,h]=CorrMatDat[i,h] * ( t(CMat[z,]-Margin[z,i]*DMat[z,])%*%M%*%(CMat[w,]-Margin[w,h]*DMat[w,]) ) /
                 sqrt( (t(CMat[z,]-Margin[z,i]*DMat[z,])%*%M%*%(CMat[z,]-Margin[z,i]*DMat[z,])) * 
                       (t(CMat[w,]-Margin[w,h]*DMat[w,])%*%M%*%(CMat[w,]-Margin[w,h]*DMat[w,])) ) }
    }
    Rrow <- cbind(Rrow,Rpart)
  }
  R <- rbind(R, Rrow)                                               # correlation matrix for test.stat
}
diag(R) <- 1


test.stat <- p.val.adj <- p.val.raw <- matrix(nrow=ncomp, ncol=nep) # matrices of test statistics and p.vals
for (z in 1:ncomp) { for (i in 1:nep) {
  test.stat[z,i]=( t(CMat[z,]-Margin[z,i]*DMat[z,])%*%meanmat[,i] ) /
                 sqrt( diag(CovMatDat)[i] * ( t(CMat[z,]-Margin[z,i]*DMat[z,])%*%M%*%(CMat[z,]-Margin[z,i]*DMat[z,]) ) )
  if (alternative=="greater") {
    p.val.adj[z,i]=1-pmvt(lower=-Inf,upper=rep(test.stat[z,i],times=ncomp*nep),df=defr,corr=R)[1]
    p.val.raw[z,i]=pt(q=test.stat[z,i],df=defr,lower.tail=FALSE) }
  if (alternative=="less") {
    p.val.adj[z,i]=1-pmvt(lower=rep(test.stat[z,i],times=ncomp*nep),upper=Inf,df=defr,corr=R)[1]
    p.val.raw[z,i]=pt(q=test.stat[z,i],df=defr,lower.tail=TRUE) }
  if (alternative=="two.sided") {
    p.val.adj[z,i]=1-pmvt(lower=rep(-abs(test.stat[z,i]),times=ncomp*nep),
                   upper=rep(abs(test.stat[z,i]),times=ncomp*nep),df=defr,corr=R)[1]
    p.val.raw[z,i]=min(pt(q=abs(test.stat[z,i]),df=defr,lower.tail=FALSE)*2,1) }
}}

list(estimate=estimate, statistic=test.stat, p.val.raw=p.val.raw, p.val.adj=p.val.adj,
     CovMatDat=CovMatDat, CorrMatDat=CorrMatDat, CorrMatComp=R, degr.fr=defr,
     Num.Contrast=CMat, Den.Contrast=DMat, Margin=Margin, alternative=alternative)


}
