'cov.LogConcDEAD' <- function (lcd) {
  if (class(lcd) != "LogConcDEAD") {
      stop("error: lcd must be of class LogConcDEAD")
  }
  triang <- lcd$triang
  x <- lcd$x
  d<-ncol(x)
  detA <- lcd$detA
  logMLE  <- lcd$logMLE
  ntriang <- nrow(triang)
  mat <- matrix(rep(0,d*d),d,d)
  for (i in 1:ntriang) {
    logMLE_triang<-logMLE[triang[i,]]
    x_triang<-matrix(x[triang[i,],],nrow=d+1)
    mat_unit_triang <- matrix(rep(0,(d+1)*(d+1)),d+1,d+1)
    for (j in 1:(d+1))
    for (k in j:(d+1)){
      mat_unit_triang[j,k]=JAD(c(logMLE_triang[c(j,k)],logMLE_triang))
    }
    mat_unit_triang = mat_unit_triang + t(mat_unit_triang)

    for (j in 1:d)
    for (k in j:d){
      coef = matrix(x_triang[,j],ncol=1) %*% matrix(x_triang[,k],nrow=1)
      mat[j,k]=mat[j,k]+sum(coef*mat_unit_triang)*detA[i]
    }
  }
 
  mat_2mom = t(mat)+mat-diag(diag(mat),nrow=nrow(mat))
  m<-matrix(colMeans(as.data.frame(x)),nrow=1)
  mat_cov<-mat_2mom- t(m)%*%m
  return (mat_cov)
}

'hatA' <- function (lcd) {
  return (cov(lcd$x)-cov.LogConcDEAD(lcd))
}












