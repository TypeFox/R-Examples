coxscor <- function(x, y, ic, offset = rep(0., length(y))) {
  ## computes cox scor function for rows of nx by n matrix  x
  ## first put everything in time order
  n <- length(y)
  nx <- nrow(x)
  yy <- y + (ic == 0.) * (1e-05)
  otag <- order(yy)
  y <- y[otag]
  ic <- ic[otag]
  x <- x[, otag, drop = F]
  ##compute  unique failure times, d=# of deaths at each failure time,
  ##dd= expanded version of d to length n, s=sum of covariates at each
  ## failure time, nn=#obs in each risk set, nno=sum(exp(offset)) at each failure time
  offset <- offset[otag]
  a <- coxstuff(x, y, ic, offset = offset)
  nf <- a$nf
  fail.times <- a$fail.times
  s <- a$s
  d <- a$d
  dd <- a$dd
  nn <- a$nn
  nno <- a$nno
  w <- rep(0., nx)
  for(i in (1.:nf)) {
    w <- w + s[, i]
    oo<- (1.:n)[y >= fail.times[i]]
    r<-rowSums(x[, oo, drop = F] * exp(offset[oo]))
    w<- w - (d[i]/nno[i])*r
  }
  return(list(scor = w, coxstuff.obj = a))
}

