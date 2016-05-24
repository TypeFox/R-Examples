probgpcm <-
function(theta,a,cb) {
    if (any(is.na(theta))) stop("theta is empty")
    if (a<=0 || is.na(a)) stop("slope is missing or negative")
    if (any(is.na(cb))) {
      ncat<-sum(!is.na(cb))+1
      if (any(is.na(cb[1:(ncat-1)]))) stop("cb is invalid")
    } else ncat<-length(cb)+1
    if (ncat<2) stop("cb is empty")
    nq<-length(theta)
    pp<-matrix(NA,nq,ncat)
    cb<-c(0,unlist(cb))
    zz<-matrix(0,nq,ncat)
    sdsum<-0
    den<-rep(0,nq)
    for (k in 1:ncat) {
      sdsum<-sdsum+cb[k]
      zz[,k]<-exp(a*(k*theta-sdsum))
      den<-den+zz[,k]
    }
    for (k in 1:ncat) {
      pp[,k]<-zz[,k]/den
    }
    return(pp)
  }
