#' @title Gumbel ES
#'
#' @description Estimates the ES of a portfolio assuming extreme losses are 
#' Gumbel distributed, for specified confidence level and holding period. 
#' Note that the long-right-hand tail is fitted to losses, not profits.
#' 
#' @param mu Location parameter for daily L/P
#' @param sigma Assumed scale parameter for daily L/P
#' @param n Assumed block size from which the maxima are drawn
#' @param cl VaR confidence level
#' @param hp VaR holding period
#' @return Estimated ES. If cl and hp are scalars, it returns scalar VaR. If cl
#' is vector and hp is a scalar, or viceversa, returns vector of VaRs. If both 
#' cl and hp are vectors, returns a matrix of VaRs.
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' National Institute of Standards and Technology, Dataplot Reference Manual. 
#' Volume 1: Commands. NIST: Washington, DC, 1997, p. 8-67.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Gumber ES Plot
#'    GumbelES(0, 1.2, 100, c(.9,.88, .85, .8), 280)
#'    
#' @export
GumbelES<- function(mu, sigma, n, cl, hp){
  
  # Check that inputs have correct dimensions
  if (!length(mu) == 1) {
    stop("Mean must be a scalar")
  }
  if (!length(sigma) == 1) {
    stop("Standard Deviation must be a scalar")
  }
  if (!is.vector(cl)) {
    stop("cl must be a vector")
  }
  if (!length(hp) == 1) {
    stop("hp must be a scalar")
  }
  
  # To check that cl is read as a scalar or row vector and hp is read as a scalar of column vector as required
  cl <- t(as.matrix(cl))
  hp <- as.matrix(hp)
  
  # Check that inputs obey sign and value restrictions
  if (sigma < 0) {
    stop("Standard deviation must be non-negative")
  }
  if (max(cl) >= 1) {
    stop("Confidence levels must be less than 1")
  }
  if (min(cl) <= 0) {
    stop("Confidence levels must be less than 1")
  }
  if (min(hp) <= 0) {
    stop("Holding period must be greated than 0")
  }
  
  # VaR Estimation
  VaR <- mu * matrix(1, dim(cl)[1], dim(cl)[2]) - sigma * log(- n * log(cl)); # Gumberl VaR
  
  # ES Estimation
  number.slices <- 1000 # Number of slices into which tail is divided
  cl0 <- cl # Initial confidence level
  term <- VaR
  
  delta.cl <- (1 - cl)/number.slices # Increment to confidence level as each slice is taken
  for (i in 1:(number.slices - 1)) {
    cl <- cl0 + i * delta.cl # Revised cl
    term <- term + mu * matrix(1, dim(cl)[1], dim(cl)[2]) - sigma * log(-n * log(cl)) # NB Gumber term
  }
  y <- term / (number.slices - 1)
  
  return(y)
}