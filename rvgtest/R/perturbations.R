## --------------------------------------------------------------------------
##
## Perturbate the given random variate generator.
##
## These generators can be used to introduce artificial errors for
## testing the test suite.
##
## --------------------------------------------------------------------------

pertadd <- function(n, rdist=rnorm, ..., min=0, max=1, p=0.001)

  ## ------------------------------------------------------------------------
  ## Function to generate random variates from a mixture of
  ## parent distribution and uniform distribution
  ## [ Add additional points uniformly on given interval (min,max) ]
  ## ------------------------------------------------------------------------
  ##  n    : The size of random sample
  ##  rdist  : Parent distribution variate generator
  ##  ...  : Parameters of the parent distribution
  ##  min  : Left boundary of uniform distribution
  ##  max  : Right boundary of uniform distribution
  ##  p    : Probability of adding error
  ## ------------------------------------------------------------------------
{ 
  x <- rdist(n,...)
  u <- runif(n,min=min,max=max)
  i <- rbinom(n,size=1,1-p)
  r <- i*x + (1-i)*u
  return(r)
}

## --------------------------------------------------------------------------

pertsub <- function(n, rdist=rnorm, ..., min=0, max=1, p=0.001)

  ## ------------------------------------------------------------------------
  ## Function to generate random variates from parent distribution
  ## rejecting points in interval (min,max) with probability p.
  ## ------------------------------------------------------------------------
  ##  n     : The size of random sample
  ##  rdist : Parent distribution variate generator
  ##  ...   : Parameters of the parent distribution
  ##  min   : Left boundary of interval
  ##  max   : right boundary of interval
  ##  p     : Probability of rejection
  ## ------------------------------------------------------------------------
{
  x <- rdist(n,...)
  i <- which(x>=min & x<=max, arr.ind=TRUE)
  v <- runif(length(i))
  j <- which(v<=p, arr.ind=TRUE)
  y <- rdist(length(j),...)
  x[i[j]] <- y
  return(x)
}

## --------------------------------------------------------------------------
