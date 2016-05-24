##
##  Difference Sign Test
##
difference.sign.test <- function(x, alternative="two.sided"){
  # Performs the difference sign test of randomness.
  #
  # Args:
  #   x: data values.
  #
  # Returns:
  #   statistic: the (normalized) value of the statistic test.
  #   n: the sample size, after the remotion of consecutive duplicate values.
  #   p.value: the asymptotic p-value.
  #
  dname <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  # Remove NAs
  x<-na.omit(x)
  n0 <- length(x)
  if (n0 < 2){stop("sample size must be greater than 1")} 
  df<-diff(x)
  # Remove differences equal to zero (from consecutive equal values)
  df<-df[df!=0]
  n <- length(df)+1
  mu <- (n-1)/2
  vr <- (n+1)/12
  test.sum <- length(df[df>0])
  test <- (test.sum-mu)/sqrt(vr)
  # p-value
  pv0 <- pnorm(test)
  if (alternative=="two.sided"){pv <- 2*min(pv0,1-pv0); alternative<-"nonrandomness"}
  if (alternative=="left.sided"){pv <- pv0; alternative<-"decreasing trend"}
  if (alternative=="right.sided") {pv <- 1-pv0; alternative<-"increasing trend"}
  # output
  #names(n)="n"
  #names(mu)="mu"
  #names(vr)="var"
  rval <- list(statistic = c(statistic=test), ds=test.sum, mu=mu, var=vr, p.value = pv, method = "Difference Sign Test", 
               data.name = dname, parameter=c(n=n), alternative=alternative) 
  class(rval) <- "htest"
  return(rval)
}

