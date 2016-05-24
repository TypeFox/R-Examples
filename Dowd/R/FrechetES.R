#' Frechet Expected Shortfall
#'
#' Estimates the ES of a portfolio  assuming extreme losses
#' are Frechet distributed, for specified confidence level and a given 
#' holding period.
#'
#' Note that the long-right-hand tail is fitted to losses, not profits.
#' 
#'
#' @param mu Location parameter for daily L/P
#' @param sigma Scale parameter for daily L/P
#' @param tail.index Tail index
#' @param n Block size from which maxima are drawn
#' @param cl Confidence level
#' @param hp Holding period
#' @return Estimated ES. If cl and hp are scalars, it returns scalar VaR. If cl
#' is vector and hp is a scalar, or viceversa, returns vector of VaRs. If both 
#' cl and hp are vectors, returns a matrix of VaRs.
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Embrechts, P., Kluppelberg, C. and Mikosch, T., Modelling Extremal Events for
#' Insurance and Finance. Springer, Berlin, 1997, p. 324.
#' 
#' Reiss, R. D. and Thomas, M. Statistical Analysis of Extreme Values from 
#' Insurance, Finance, Hydrology and Other Fields, Birkhaueser, Basel, 1997, 
#' 15-18.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes ES assuming Frechet Distribution for given parameters
#'    FrechetES(3.5, 2.3, 1.6, 10, .95, 30)
#'
#' @export
FrechetES <- function(mu, sigma, tail.index, n, cl, hp){
  
  # Check that inputs have correct dimensions
  if (!length(mu) == 1) {
    stop("mu must be a scalar")
  }
  if (!length(sigma) == 1) {
    stop("sigma must be a scalar")
  }
  if (!length(tail.index) == 1) {
    stop("tail.index must be a scalar")
  }
  if (!is.vector(cl)) {
    stop("cl must be a vector or a scalar")
  }
  if (!is.vector(hp)) {
    stop("hp must be a vector or a scalar")
  }
  
  # Change cl and hp to row vector and column vectors respectively
  cl <- t(as.matrix(cl))
  hp <- as.matrix(hp)
  
  # Check that parameters obey sign and value restrictions
  if (sigma < 0) {
    stop("Standard deviation must be non-negative")
  }
  if (min(tail.index) <= 0) {
    stop("Tail index must be greater than 0")
  }
  if ( max(cl) >= 1){
    stop("Confidence level(s) must be less than 1")
  }
  if ( min(cl) <= 0){
    stop("Confidence level(s) must be greater than 0")
  }
  if ( min(cl) <= 0){
    stop("Holding period(s) must be greater than 0")
  }
  
  # VaR estimation
  VaR <- mu * matrix(1, 1, length(cl)) - (sigma / tail.index) * 
    (1 - ( - n * log(cl)) ^ ( - tail.index))
  
  # ES Estimation
  number.slices <- 1000 # Number of slices into which tail is divided
  cl0 <- cl # Initial confidence level
  term <- VaR
  
  delta.cl <- (1 - cl) / number.slices # Increment to confidence level as each slice is taken
  for (i in 1:(number.slices-1)) {
    cl <- cl0 + i * delta.cl # Revised cl
    term <- term + mu * matrix(1, 1, length(cl)) - (sigma / tail.index) * 
      (1 - ( - n * log(cl)) ^ ( - tail.index))
    # NB Frechet term
  }
  
  y <- term / (number.slices - 1)
  
  return(y)
  
} 