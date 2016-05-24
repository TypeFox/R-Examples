##' Does a simple test for the median based on the binomial distribution
##'
##' @param x data
##' @param median value to test
##' @return p-value
##' @export
simple.median.test <- function (x,median=NA) {
  ### does a simle test for the median based on the binomial distribution
  ### of the number of observations larger (or smaller) than the true median.
  x <- as.vector(x)
  n <- length(x)
  ## hypothesis test, two-sided
  bigger<-sum(x>median)
  equal <- sum(x==median)
  count <- bigger + equal/2
  count <- min(count,n-count)
  ret <- 2*pbinom(count,n,.5)

  ret
}
