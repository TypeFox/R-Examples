
coxvar <- function(x, y, ic, offset = rep(0., length(y)), coxstuff.obj = NULL){
  ## computes information elements (var) for cox
  ## x is nx by n matrix of expression  values
  nx <- nrow(x)
  n <- length(y)
  yy <- y + (ic == 0.) * (1e-06)
  otag <- order(yy)
  y <- y[otag]
  ic <- ic[otag]
  x <- x[, otag, drop = F]
  offset <- offset[otag]
  if(is.null(coxstuff.obj)) {
    coxstuff.obj <- coxstuff(x, y, ic, offset = offset)
  }
  nf <- coxstuff.obj$nf
  fail.times <- coxstuff.obj$fail.times
  s <- coxstuff.obj$s
  d <- coxstuff.obj$d
  dd <- coxstuff.obj$dd
  nn <- coxstuff.obj$nn
  nno <- coxstuff.obj$nno
  x2<- x^2
  oo <- (1.:n)[y >= fail.times[1] ]
  sx<-(1/nno[1])*rowSums(x[, oo] * exp(offset[oo]))
  s<-(1/nno[1])*rowSums(x2[, oo] * exp(offset[oo]))
  w <-  d[1] * (s - sx * sx)
  for(i in 2.:nf) {
    oo <- (1.:n)[y >= fail.times[i-1] & y < fail.times[i] ]
    sx<-(1/nno[i])*(nno[i-1]*sx-rowSums(x[, oo,drop=F] * exp(offset[oo])))
    s<-(1/nno[i])*(nno[i-1]*s-rowSums(x2[, oo,drop=F] * exp(offset[oo])))
    w <- w + d[i] * (s - sx * sx)
  }
  return(w)
}
