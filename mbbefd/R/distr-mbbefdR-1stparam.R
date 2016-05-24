

### R version of d,p,q,r functions MBBEFD(a,b)

dmbbefdR <- function(x, a, b, log=FALSE, g)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(x)))
  if(!missing(g)) #to be removed
    stop("please use dMBBEFD for the second parametrization.")
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- rep(0, length(x))
    res[x == 1] <- 1
  }else if(is.infinite(a))
  {
    res <- -log(b)*b^x
  }else
  {
    res <- -a * (a+1) * b^x * log(b) / (a + b^x)^2 
    res[x == 1] <- (a+1) * b / (a+b)
  }
  res[x > 1] <- 0
  res[x < 0] <- 0
  
  if(log)
    res <- log(res)
  res  
}  
	
pmbbefdR <- function(q, a, b, lower.tail = TRUE, log.p = FALSE, g)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(q)))
  if(!missing(g)) #to be removed
    stop("please use pMBBEFD for the second parametrization.")
  
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- rep(0, length(q))
  }else if(is.infinite(a))
  {
    res <- 1-b^q
  }else
  {
    res <- a * ( (a+1) / (a + b^q) - 1) 
  }
  res[q >= 1] <- 1
  res[q < 0] <- 0

  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  
  res
}  
  
	 
qmbbefdR <- function(p, a, b, lower.tail = TRUE, log.p = FALSE, g)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(p)))
  if(!missing(g)) #to be removed
    stop("please use qMBBEFD for the second parametrization.")
  
  
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p) 
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- rep(0, length(p))
    res[p > 0] <- 1
  }else if(is.infinite(a))
  {
    res <- rep(1, length(p))
    res[p < 1-b] <- log(1-p[p < 1-b])/log(b)
  }else
  {
    pab <- (a+1)*b/(a+b)
    res <- rep(1, length(p))
    p2 <- p[p < 1-pab]
    res[p < 1-pab] <- log((1-p2)*a/(a+p2))/log(b)
  }
  res[p < 0 | p > 1] <- NaN
  
  res
}  

  
rmbbefdR <- function(n, a, b)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, n))
  qmbbefdR(runif(n, 0, 1), a, b)
}
	
	
ecmbbefdR <- function(x, a, b)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(x)))
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- x
  }else if(is.infinite(a))
  {
    res <- (1-b^x)/(1-b)
  }else
  {
    res <- log((a+b^x)/(a+1))/log((a+b)/(a+1))  
  }
  res[x < 0] <- 0
  res[x > 1] <- 1
  
  res
}
	
#moment
mmbbefdR <- function(order, a, b)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(order)))
  
  if(order == 1)
    return(log((a+b)/(a+1))/log(b)*(a+1))
  else if(order == 2)
  {
    2*(a+1)/log(b)*(log(a+b) - gendilog(a,b))
  }else
    stop("not yet implemented.")
}
	
#total loss
tlmbbefdR <- function(a, b)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(NaN)
  if(is.infinite(a))
    res <- b
  else
    res <- (a+1)*b/(a+b)
  res
}



### d,p,q,ec,m,tl functions MBBEFD(a,b)

dmbbefd <- dmbbefdR

pmbbefd <- pmbbefdR

qmbbefd <- qmbbefdR

ecmbbefd <- ecmbbefdR

mmbbefd <- mmbbefdR

tlmbbefd <- tlmbbefdR

