
##  Part of R package BDA
##  Copyright (C) 2009-2010 Bin Wang
##
##  Unlimited use and distribution (see LICENCE).

##  2012/10/31: define four functions for a mixture of k normal
##  components f(x) = \sum p_i * dnorm(x,mu_i, sigma_i)

.dmnorm <- function(x,p,mean,sd){
  k <- length(p)
  res <- 0
  for(i in 1:k){
    res <- res + p[i] * dnorm(x,mean[i],sd[i])
  }
  res
}

dmixnorm <- function(x,p,mean,sd){
  if(missing(p)) p <- 1
  if(missing(mean)) mean <- 0
  if(missing(sd)) sd <- 1
  ndim <- length(p)
  if(length(mean) != ndim | length(sd) != ndim)
    stop("Parameters have different lengths")
  if(any(p>1|p<0)) stop("Wrong mixing coefficients")
  p <- p/sum(p)
  if(any(sd<=0)) stop("Invalid standard deviation(sd)")
  
  sapply(x,.dmnorm,p=p,mean=mean,sd=sd)
}

.pmnorm <- function(x,p,mean,sd){
  k <- length(p)
  res <- 0
  for(i in 1:k){
    res <- res + p[i] * pnorm(x,mean[i],sd[i])
  }
  res
}

pmixnorm <- function(q,p,mean,sd){
  if(missing(p)) p <- 1
  if(missing(mean)) mean <- 0
  if(missing(sd)) sd <- 1
  ndim <- length(p)
  if(length(mean) != ndim | length(sd) != ndim)
    stop("Parameters have different lengths")
  if(any(p>1|p<0)) stop("Wrong mixing coefficients")
  p <- p/sum(p)
  if(any(sd<=0)) stop("Invalid standard deviation(sd)")
  
  sapply(q,.pmnorm,p=p,mean=mean,sd=sd)
}

.rmnorm <- function(x,p){
  k <- length(p)
  cump <- cumsum(p)
  x[which(runif(1)-cump<=0)[1]]
}

rmixnorm <- function(n,p,mean,sd){
  if(missing(p)) p <- 1
  if(missing(mean)) mean <- 0
  if(missing(sd)) sd <- 1
  ndim <- length(p)
  if(length(mean) != ndim | length(sd) != ndim)
    stop("Parameters have different lengths")
  if(any(p>1|p<0)) stop("Wrong mixing coefficients")
  p <- p/sum(p)
  if(any(sd<=0)) stop("Invalid standard deviation(sd)")
  n <- ceiling(n)
  stopifnot(n>0)
  
  tmp <- NULL
  k <- length(p)
  for(i in 1:k){
    tmp <- cbind(tmp, rnorm(n,mean[i], sd[i]))
  }
  res <- apply(tmp,1,.rmnorm,p=p)
  as.numeric(res)
}

qmixnorm <- function(prob,p,mean,sd){
  if(missing(p)) p <- 1
  if(missing(mean)) mean <- 0
  if(missing(sd)) sd <- 1
  sele <- !is.na(prob)
  if(any(prob[sele]>1|prob[sele]<0))
    stop("Invalid 'prob' value(s)")

  ndim <- length(p)
  if(length(mean) != ndim | length(sd) != ndim)
    stop("Parameters have different lengths")
  if(any(p>1|p<0)) stop("Wrong mixing coefficients")
  p <- p/sum(p)
  if(any(sd<=0)) stop("Invalid standard deviation(sd)")

  mean.pool <- sum(p*mean)
  s.pool <- sqrt(sum(p^2*sd^2))  
  x <- seq(mean.pool-4*s.pool,mean.pool+4*s.pool,length=401L)
  Fx <- sapply(x,.pmnorm,p=p,mean=mean,sd=sd)
  approx(Fx,x,prob)$y
}
