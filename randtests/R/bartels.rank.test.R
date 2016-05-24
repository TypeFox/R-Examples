##
##  Bartels' Rank Test
##
bartels.rank.test <- function(x, alternative="two.sided", pvalue="normal"){
  # Performs Bartels Ratio Test for Randomness.
  #
  # Args:
  #   x: a numeric vector containing the data.
  #   alternative hypothesis, must be one of "two.sided" (default), "left.sided" or "right.sided"
  #   pv.method: asymptotic aproximation method used to compute the p-value.
  #
  # Returns:
  #   statistic: the value of the RVN statistic test and the theoretical mean value and variance of the RVN statistic test.
  #   n: the sample size, after the remotion of consecutive duplicate values.
  #   p.value: the asymptotic p-value.
  #   method: a character string indicating the test performed.
  #   data.name: a character string giving the name of the data.
  #   alternative: a character string describing the alternative.  
  #
  dname <- deparse(substitute(x))
  # Remove NAs
  x<-na.omit(x)
  stopifnot(is.numeric(x))
  n <- length(x)
  if (alternative == "t"){alternative <- "two.sided"} 
  if (alternative == "l"){alternative <- "left.sided"}
  if (alternative == "r"){alternative <- "right.sided"}    
  if (alternative != "two.sided" & alternative != "left.sided" & alternative != "right.sided")
    {stop("must give a valid alternative")}
  if (n < 3){stop("sample size must be greater than 2")}
  # unique
  rk <- rank(x)
  d <- diff(rk)
  #d.rank <- n*(n^2-1)/12
  nm <- sum(rk^2)
  d.rank <- nm-n*(mean(rk)^2)
  RVN <- sum(d^2)/d.rank
  mu <- 2
  vr <- (4*(n-2)*(5*n^2-2*n-9))/(5*n*(n+1)*(n-1)^2)
  # 
  # Computes the p-value
  if (pvalue == "auto"){
    pvalue <- ifelse(n<=100,"beta","normal")
    if (n<10) pvalue <- "exact"
  }
  if (pvalue == "exact"){
    pv <- pbartelsrank(q=c(nm-1,nm), n)
    pv0 <- pv[2]
    pv1 <- pv[1]
  }
    
  if (pvalue == "beta"){
     btp <- (5*n*(n+1)*(n-1)^2)/(2*(n-2)*(5*n^2-2*n-9))-1/2
     pv0 <- pbeta(RVN/4,shape1=btp,shape2=btp)
     pv1 <- 1-pv0
  }
  if (pvalue=="normal") pv0 <- pnorm((RVN - mu) / sqrt(vr))
  if (pvalue=="beta" | pvalue=="normal") pv1 <- 1-pv0  
  #
  if (alternative=="two.sided"){
    pv <- 2*min(pv0, pv1) 
    alternative<-"nonrandomness"
  }
  if (alternative=="left.sided"){
    pv <- pv0 
    alternative<-"trend"
  }
  if (alternative=="right.sided"){
    pv <- pv1
    alternative<-"systematic oscillation"
  }

  test <- (RVN - mu) / sqrt(vr)
  rval <- list(statistic = c(statistic=test), nm=sum(d^2), rvn=RVN, mu=mu, var=vr, p.value = pv, 
               method = "Bartels Ratio Test", data.name = dname, parameter=c(n=n), n=n, alternative=alternative)  
  class(rval) <- "htest"
  return(rval)
  
}
