
#<<BEGIN>>
dbetagen <- function(x,shape1,shape2,min=0,max=1,ncp=0,log=FALSE)
#TITLE The Generalised Beta Distribution
#NAME betagen
#KEYWORDS distribution
#DESCRIPTION
#Density, distribution function, quantile function and random generation for the Beta distribution
#defined on the \samp{[min, max]} domain with parameters \samp{shape1} and \samp{shape2} (
#and optional non-centrality parameter \samp{ncp}).
#INPUTS
#{x,q}<<Vector of quantiles.>>
#{p}<<Vector of probabilities.>>
#{n}<<Number of observations. If \samp{length(n) > 1}, the length is taken to be the number required.>>
#{shape1, shape2}<<Positive parameters of the Beta distribution.>>
#[INPUTS]
#{min}<<Vector of minima.>>
#{max}<<Vector of maxima.>>
#{ncp}<<Non-centrality parameter of the Beta distribution.>>
#{log, log.p}<<Logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.>>
#{lower.tail}<<Logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.>>
#DETAILS
#\deqn{x \sim  betagen(shape1, shape2, min, max, ncp)}{x ~ betagen(shape1, shape2, min, max, ncp)}
#if
#\deqn{\frac{x-min}{max-min}\sim beta(shape1,shape2,ncp)}{(x-min)/(max-min)~beta(shape1,shape2,ncp)}
#These functions use the \code{\link{Beta}} distribution functions after correct parametrisation.
#EXAMPLE
#curve(dbetagen(x,shape1=3,shape2=5,min=1,max=6), from = 0, to = 7)
#curve(dbetagen(x,shape1=1,shape2=1,min=2,max=5), from = 0, to = 7, lty=2, add=TRUE)
#curve(dbetagen(x,shape1=.5,shape2=.5,min=0,max=7), from = 0, to = 7, lty=3, add=TRUE)
#SEE ALSO
#\code{\link{Beta}}
#VALUE
#\samp{dbetagen} gives the density, \samp{pbetagen} gives the distribution function,
#\samp{qbetagen} gives the quantile function, and \samp{rbetagen} generates random deviates.

#CREATED 08-04-16
#--------------------------------------------
{
  if(length(x) == 0) return(numeric(0))
  ow <- options(warn=-1)
    min <- as.vector(min)
	max <- as.vector(max)

  x <- (x - min)/(max - min)
  options(ow)
  if(missing(ncp))
  d <- dbeta(x, shape1=shape1, shape2=shape2) / (max-min)
  else
  {
  ncp<- as.vector(ncp)
  d <- dbeta(x, shape1=shape1, shape2=shape2, ncp=ncp) / (max-min)
  }
  if(log) d <- log(d)
  d[max <= min] <- NaN
  if(any(is.na(d))) warning("NaN in dbetagen")
  return(d)}

#<<BEGIN>>
pbetagen <- function(q,shape1,shape2,min=0,max=1,ncp=0,lower.tail = TRUE, log.p = FALSE)
#ISALIAS dbetagen
#--------------------------------------------
{
  if(length(q) == 0) return(numeric(0))
    min <- as.vector(min)
	max <- as.vector(max)
  q2 <- (q - min)/(max-min)
  ow <- options(warn=-1)
  if(missing(ncp))
  p <- pbeta(q2,shape1=shape1,shape2=shape2, lower.tail=lower.tail,log.p=log.p)
  else
  {
  ncp <- as.vector(ncp)
  p <- pbeta(q2,shape1=shape1,shape2=shape2, ncp=ncp, lower.tail=lower.tail, log.p=log.p)
  }
  options(ow)
  # If min = max = q -> should return 1
  quel <- (abs(q - min) < (.Machine$double.eps^0.5)) & 
          (abs(q - max) < (.Machine$double.eps^0.5))  #if min == max == q
  p[quel] <- if(lower.tail) 1 else 0 
  if(log.p) p[quel] <- log(p[quel]) 
    
  p[max < min] <- NaN
  if(any(is.na(p))) warning("NaN in pbetagen")
  return(p)}

#<<BEGIN>>
qbetagen <- function(p,shape1,shape2,min=0,max=1,ncp=0,lower.tail=TRUE,log.p=FALSE)
#ISALIAS dbetagen
#--------------------------------------------
{
  if(length(p) == 0) return(numeric(0))
    min <- as.vector(min)
	max <- as.vector(max)
  ow <- options(warn=-1)
  if(missing(ncp))
  q <- qbeta(p,shape1=shape1, shape2=shape2, lower.tail=lower.tail, log.p=log.p)
  else
  {
  ncp <- as.vector(ncp)
  q <- qbeta(p,shape1=shape1, shape2=shape2, ncp=ncp, lower.tail=lower.tail, log.p=log.p)
  }
  options(ow)
  q2 <- q * (max-min) + min
  q2[max < min] <- NaN
  if(any(is.na(q2))) warning("NaN in qbetagen")
  return(q2)}


#<<BEGIN>>
rbetagen <- function(n,shape1,shape2,min=0,max=1,ncp=0)
#ISALIAS dbetagen
#--------------------------------------------
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rbetagen")  
    min <- as.vector(min)
	max <- as.vector(max)
  ow <- options(warn=-1)
  if(missing(ncp))
  r <- rbeta(n, shape1=shape1, shape2=shape2)
  else
  {
  ncp <- as.vector(ncp)
  r <- rbeta(n, shape1=shape1, shape2=shape2, ncp=ncp)
  }
  options(ow)
  r2 <- r*(max-min) + min
  r2[max < min] <- NaN
  if(any(is.na(r2))) warning("NaN in rbetagen")
  return(r2)
}
