#d, p, q, r function for one-inflated distribution


doifun <- function(x, dfun, p1, log=FALSE, ...)
{
  if(!(p1 >= 0 && p1 <= 1))
    return(rep(NaN, length(x)))
  
  res <- dfun(x, log=FALSE, ...)*(1 - p1)*(x != 1) + p1*(x == 1)
  
  if(log)
    res <- log(res) 
  res
}

poifun <- function(q, pfun, p1, lower.tail = TRUE, log.p = FALSE, ...)
{
  if(!(p1 >= 0 && p1 <= 1))
    return(rep(NaN, length(q)))
  
  res <- pfun(q, lower.tail = TRUE, log.p = FALSE, ...)*(1 - p1) + p1*(q >= 1)
  
  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  
  res
}


qoifun <- function(p, qfun, p1, lower.tail = TRUE, log.p = FALSE, ...)
{
  if(!(p1 >= 0 && p1 <= 1))
    return(rep(NaN, length(p)))
  
  p <- p/(1-p1) #transformed quantile
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p) 
  
  res <- qfun(p, lower.tail = TRUE, log.p = FALSE, ...)
  res[p >= 1-p1] <- 1
  
  res
}  

roifun <- function(n, rfun, p1, ...)
{
  n <- ifelse(length(n)>1, length(n), n)
  if(!(p1 >= 0 && p1 <= 1))
    return(rep(NaN, n))
  res <- rfun(n, ...)
  res[rbinom(n, 1, p1) == 1] <- 1
  res
}

#exposure curve and moment functions
ecoifun <- function(x, ecfun, mfun, p1, ...)
{
  if(!(p1 >= 0 && p1 <= 1))
    return(rep(NaN, length(x)))
  
  G0 <- ecfun(x, ...) #exposure curve
  E0 <- mfun(order=1, ...) #expectation
  
  ((1-p1)*G0 + p1*x/E0)/(1-p1+p1/E0)
}


# moment function
moifun <- function(order, mfun, p1, ...)
{
  if(!(p1 >= 0 && p1 <= 1))
    return(rep(NaN, length(order)))
  
  E0 <- mfun(order=order, ...) #expectation
  p1 + (1-p1)*E0
}

#total loss function
tloifun <- function(p1, ...)
{
  p1
}
