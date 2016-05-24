##
##  Mann-Kendall Rank Test
##
rank.test <- function(x, alternative="two.sided")
{
  dname <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  # Remove NAs
  x<-na.omit(x)
  
  # Main code  
  n <- length(x)
  p <- 0
  for (i in 1:(n-1)){
    t <- x[(i+1):n]
    p <- p+length(t[t>x[i]])
  }
  mu <- n*(n-1)/4
  vr <- n*(n-1)*(2*n+5)/72
  test <- (p-mu)/sqrt(vr)
  # p-value
  pv0 <- pnorm(test)
  if (alternative=="two.sided"){
    pv <- 2*min(pv0,1-pv0) 
    alternative<-"trend"
  }
  if (alternative=="left.sided"){
    pv <- pv0 
    alternative<-"downnward trend"
  }
  if (alternative=="right.sided"){
    pv <- 1-pv0 
    alternative<-"upward trend"
  }  
  # Output
  rval <- list(statistic = c(statistic=test), P=p, mu=mu, var=vr, p.value = pv, method = "Mann-Kendall Rank Test", 
               data.name = dname, parameter=c(n=n), alternative=alternative) 
  class(rval) <- "htest"
  return(rval)
}
