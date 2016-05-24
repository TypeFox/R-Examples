#' Pickands Estimator
#'
#' Estimates the Value of Pickands Estimator for a specified data set
#' and chosen tail size. Notes: (1) We estimate the Pickands Estimator by 
#' looking at the upper tail. (2) The tail size must be less than one quarter 
#' of the total sample size. (3) The tail size must be a scalar.
#'
#' @param Ra A data set
#' @param tail.size Number of observations to be used to estimate the Pickands
#' estimator
#' @return Value of Pickands estimator
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes estimated Pickands estimator for randomly generated data.
#'    Ra <- rnorm(1000)
#'    PickandsEstimator(Ra, 40)
#'
#' @export
PickandsEstimator <- function (Ra, tail.size) {
  
  data <- as.vector(Ra)
  data <- sort(data)
  n <- length(data)
  k <- tail.size
  num <- data[n - k] - data[n - 2 * k]
  den <- data[n - 2 * k] - data[n - 4 * k]
  y <- log(num / den) / log(2)
  return(y)
  
}