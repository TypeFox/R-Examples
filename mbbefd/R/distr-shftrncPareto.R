#d, p, q, r function for shifted truncated Pareto distribution

#should it be dstpareto01?
dstpareto <- function(x, a, log=FALSE)
{
  if(!(a > 0))
    return(rep(NaN, length(x)))
  
  res <- a * (x+1)^(-a-1)/(1 - 2^(-a))
  res[x > 1] <- 0
  res[x < 0] <- 0
  
  if(log)
    res <- log(res) 
  res
}

pstpareto <- function(q, a, lower.tail = TRUE, log.p = FALSE)
{
  if(!(a > 0))
    return(rep(NaN, length(q)))
  
  res <- (1 - (q+1)^(-a))/(1-2^(-a))
  res[q >= 1] <- 1
  res[q < 0] <- 0
  
  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  
  res
}


qstpareto <- function(p, a, lower.tail = TRUE, log.p = FALSE)
{
  if(!(a > 0))
    return(rep(NaN, length(p)))
  
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p) 
  
  res <- (1-p*(1-2^(-a)))^(-1/a) - 1
  res[p < 0 | p > 1] <- NaN
  
  res
}  

rstpareto <- function(n, a)
{
  n <- ifelse(length(n)>1, length(n), n)
  if(!(a > 0))
    return(rep(NaN, n))
  qstpareto(runif(n, 0, 1), a)
}


ecstpareto <- function(x, a)
{
  if(!(a > 0))
    return(rep(NaN, length(x)))
  
  if(a == 1)
  {
    res <- (2*log(x+1) - x)/(2*log(2) - 1)
  }else
  {
    res <- ((x+1)^(-a+1) - 2^(-a)*x*(-a+1) - 1)/(2^(-a+1)-2^(-a)*(-a+1) - 1)
  }
  res[x < 0] <- 0
  res[x > 1] <- 1
  
  res
}

mstpareto <- function(order, a)
{
  if(order == 1)
    return(ifelse(a == 1, 2*log(2)-1, (2^(-a+1) - 2^(-a)*(-a+1)-1)/(-a+1)/(1-2^(-a)) ))
  else
    stop("not yet implemented")
}
  
