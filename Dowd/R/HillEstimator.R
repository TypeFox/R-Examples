#' Hill Estimator
#'
#' Estimates the value of the Hill Estimator for a given specified data set and
#' chosen tail size. Notes: 1) We estimate Hill Estimator by looking at the 
#' upper tail. 2) If the specified tail size is such that any included 
#' observations are negative, the tail is truncated at the point before 
#' observations become negative. 3) The tail size must be a scalar.
#' 
#' @param Ra Data set
#' @param tail.size Number of observations to be used to estimate the Hill
#' estimator.
#' @return Estimated value of Hill Estimator
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates Hill Estimator of 
#'    Ra <- rnorm(15)
#'    HillEstimator(Ra, 10)
#'
#' @export
HillEstimator <- function(Ra, tail.size){
  
  data <- as.vector(Ra)
  n <- length(data)
  i <- which(data <= 0)
  i <- max(i)
  k <- min(tail.size, n - i)
  logsum <- 0
  for (i in 2:k){
    logsum <- logsum + log(data[n - 1])
  }
  y <- logsum / (k - 1) - log(data[n - k - 1])
  return(y)
  
}