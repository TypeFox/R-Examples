coxstuff<- function(x, y, ic, offset = rep(0., length(y))) {
  fail.times <- unique(y[ic == 1.])
  nf <- length(fail.times)
  n <- length(y)
  nn <- rep(0., nf)
  nno <- rep(0., nf)
  for(i in 1.:nf) {
    nn[i] <- sum(y >= fail.times[i])
    nno[i] <- sum(exp(offset)[y >= fail.times[i]])
  }
  s <- matrix(0., ncol = nf, nrow = nrow(x))
  d <- rep(0., nf)
  ##expand d out to a vector of length n
  for(i in 1.:nf) {
    o <- (1.:n)[(y == fail.times[i]) & (ic == 1.)]
    d[i] <- length(o)
  }
  oo <- match(y, fail.times)
  oo[ic==0]<-NA
  oo[is.na(oo)]<- max(oo[!is.na(oo)])+1
  s<-t(rowsum(t(x),oo))
  if(ncol(s)> nf){s<-s[,-ncol(s)]}
  dd <- rep(0., n)
  for(j in 1.:nf) {
    dd[(y == fail.times[j]) & (ic == 1.)] <- d[j]
  }
  return(list(fail.times=fail.times, s=s, d=d, dd=dd, nf=nf, nn=nn, nno=nno))
}
