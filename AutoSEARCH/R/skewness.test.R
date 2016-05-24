skewness.test <-
function(x)
{
  x.n <- length(x)
  x.mean <- mean(x)
  x.demeaned <- x - x.mean
  numer <- I(x.n^2)*sum(x.demeaned^3)^2
  denom <- 6*sum(x.demeaned^2)^3
  statistic <- numer/denom

  out <- list()
  out$statistic <- statistic
  out$p.value <- pchisq(statistic, 1, lower.tail=FALSE)

  return(out)
}
