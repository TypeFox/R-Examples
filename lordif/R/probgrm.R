probgrm <-
function(theta,a,cb) {
    if (any(is.na(theta))) stop("theta is empty")
    if (a<=0 || is.na(a)) stop("slope is missing or negative")
    if (any(is.na(cb))) {
      ncat<-sum(!is.na(cb))+1
      if (any(is.na(cb[1:(ncat-1)]))) stop("cb is invalid")
    } else ncat<-length(cb)+1
    if (ncat<2) stop("cb is empty")
    if (all(order(cb[!is.na(cb)])!=1:(ncat-1))) stop("cb is disordinal")
    nq<-length(theta)
    pp<-matrix(NA,nq,ncat)
    ps<-matrix(0,nq,ncat)
    ps[,1]<-1
    for (k in 1:(ncat-1)) {
      ps[,k+1]<-1/(1+exp(-a*(theta-unlist(cb[k]))))
    }
    pp[,1]<-1-ps[,1]
    pp[,ncat]<-ps[,ncat]
    for (k in 1:(ncat-1)) {
      pp[,k]=ps[,k]-ps[,k+1]
    }
    return(pp)
  }
