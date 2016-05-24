xBalance.makepooledsd <- function(zz,mm,pre.n) {
  if (any(zz>0)) {
    s2.t <- apply(mm[(zz>0),,drop=FALSE],2,var,na.rm=TRUE)
  } else {
    s2.t <- 0
  }

  ##Variance of the covariates among the controls (zz==0)
  if (any(zz<=0)) {
    s2.c <- apply(as.matrix(mm[(zz<=0),,drop=FALSE]), 2, var,na.rm=TRUE)
  } else {
    s2.c <- 0
  }

  ##Pooled standard deviation.
  sqrt( (max(sum(zz>0)-1,0)*s2.t + max(sum(zz<=0)-1,0)*s2.c)/(pre.n-2) )
}
