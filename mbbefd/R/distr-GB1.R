#d, p, q, r function for generalized beta distribution of the first kind
#(no location and scale paramater)

#should it be dstpareto01?
dgbeta <- function(x, shape0, shape1, shape2, log=FALSE)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(x)))
  if(log)
  {
    res <- log(shape0)+(shape0-1)*log(x)+dbeta(x^(shape0), shape1, shape2, log=log)
    res[x < 0 | x > 1] <- -Inf
  }else
  {
    res <- shape0*x^(shape0-1)*dbeta(x^(shape0), shape1, shape2, log=FALSE)
    res[x < 0 | x > 1] <- 0
  }
  res
}

pgbeta <- function(q, shape0, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(q)))
  res <- pbeta(q^shape0, shape1, shape2, lower.tail=TRUE, log.p=FALSE)
  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  res
}


qgbeta <- function(p, shape0, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(p)))
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p)
  qbeta(p, shape1, shape2, lower.tail=TRUE, log.p=FALSE)^(1/shape0)
}  

rgbeta <- function(n, shape0, shape1, shape2)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  n <- ifelse(length(n)>1, length(n), n)
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, n))
  rbeta(n, shape1, shape2)^(1/shape0)
}


ecgbeta <- function(x, shape0, shape1, shape2)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(x)))
  
  cst2 <- beta(shape1, 1/shape0)/beta(shape1 + shape2, 1/shape0)
  
  pbeta(x^shape0, shape1+1/shape0, shape2) + x*(1 - pbeta(x^shape0, shape1, shape2))*cst2
}

mgbeta <- function(order, shape0, shape1, shape2)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(order)))
  
  beta(shape1+shape2, order/shape0) / beta(shape1, order/shape0)
}
  

###################
#internal functions

#incomplete beta function
betainc <- function(x, a,b) pbeta(x, a, b)*beta(a,b)


#Theil index, see package ineq for other income index (e.g. Gini coefficient)
Theil.theo  <- function(shape0, shape1, shape2)
{
  EX <- beta(shape1+shape2, 1/shape0) / beta(shape1, 1/shape0)
  1/shape0*(digamma(shape1+1/shape0)-digamma(shape1+shape2+ 1/shape0)) - log(EX)
}

Theil.theo.shape0  <- function(shape0, obs)
{
  #compute shape1/shape2 on a rescaled sample and moment estimator
  obs <- obs^shape0
  n <- length(obs)
  m <- mean(obs)
  v <- (n - 1)/n*var(obs)
  aux <- m*(1-m)/v - 1
  shape1 <- m*aux
  shape2 <- (1-m)*aux
  
  Theil.theo(shape0, shape1, shape2)
}
