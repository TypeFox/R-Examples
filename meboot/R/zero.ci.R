
zero.ci <- function(x, confl=0.05)
{
  bigj <- length(x)
  yle <- x[x<=0]  #le=less than or equal to 
  m1 <- length(yle)
  #m1/bigj
  ygt <- x[x>0]  #gt=greater than
  m2 <- length(ygt)
  #m2/bigj

  # Want confidence interval adjusted so it covers 
  # the true zero (1-confl)*100 times one realization gives 
  # on random interval
  xsort <- sort(x)
  nlo <- ceiling(confl*m1)
  nup <- floor(confl*m2)
  if(nlo <= 0)    nlo <- 1
  if(nup == bigj) nup <- 0
  lolim <- xsort[nlo]
  uplim <- xsort[bigj-nup]
  bnlo <- length(which(x<=lolim))
  bnup <- length(which(x>uplim))

  list(bnlo=bnlo, bnup=bnup, lolim=lolim, uplim=uplim)
}
