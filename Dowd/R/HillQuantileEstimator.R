#' Hill Quantile Estimator
#'
#' Estimates value of Hill Quantile Estimator for a specified data set, tail
#' index, in-sample probability and confidence level.
#'
#' @param Ra A data set
#' @param tail.index Assumed tail index
#' @param in.sample.prob In-sample probability (used as basis for projection)
#' @param cl Confidence level
#' @return Value of Hill Quantile Estimator
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Next reference
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes estimates value of hill estimator for a specified data set
#'    Ra <- rnorm(1000)
#'    HillQuantileEstimator(Ra, 40, .5, .9)
#'
#' @export
HillQuantileEstimator <- function(Ra, tail.index, in.sample.prob, cl){
  
  data <- as.vector(Ra)
  data <- sort(data)
  n <- length(data)
  a <- in.sample.prob * n
  k <- ((a>=0)*floor(a)+(a<0)*ceiling(a))
  p <- 1 - cl
  y <- data[n - k] * (p * n / k) ^ (- tail.index)
  return(y)
  
} 