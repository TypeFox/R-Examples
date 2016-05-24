nbfactors <-
function(dta,maxnbfactors=12,diagnostic.plot=FALSE,minerr=1e-03, maxiter=100){
  dig = 2
  m = ncol(dta)
  n = nrow(dta)
  falist = vector(length=maxnbfactors+1,"list")
  falist[[1]] = list(B=matrix(0,ncol=1,nrow=m))
  falist[-1] = lapply(1:maxnbfactors,emfahighdim,eps=dta,minerr=minerr, maxiter=maxiter)
  Blist = lapply(falist,function(fa,m) matrix(fa$B,nrow=m),m=m)
  sdt = VarInflation(dta,Blist,maxnbfactors,dig)
  if (diagnostic.plot) {
    plot(0:maxnbfactors,sdt,ylab="Variance Inflation Criterion",xlab="Number of factors",bty="l",
         lwd=1.25,type="b",pch=16,cex.lab=1.25,cex=1.25,cex.axis=1.25)
  }
  opt <- which.min(sdt)-1
  list(criterion=sdt,optimalnbfactors=opt)
}
