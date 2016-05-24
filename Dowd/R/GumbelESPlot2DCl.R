#' @title Gumbel VaR
#'
#' @description Estimates the EV VaR of a portfolio assuming extreme losses are Gumbel distributed, for specified confidence level and holding period.
#'
#' @param mu Location parameter for daily L/P
#' @param sigma Assumed scale parameter for daily L/P
#' @param n size from which the maxima are drawn
#' @param cl VaR confidence level
#' @param hp VaR holding period
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots ES against Cl
#'    GumbelESPlot2DCl(0, 1.2, 100, seq(0.8,0.99,0.02), 280)
#'
#' @export
GumbelESPlot2DCl<- function(mu, sigma, n, cl, hp){
  
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
  
  # To check that cl is read as a row vector as required
  cl <- t(as.matrix(cl))
  
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
  VaR <- mu * matrix(1, dim(cl)[1], dim(cl)[2]) - sigma * log(- n * log(cl));
  
  # ES Estimation
  number.slices <- 1000 # Number of slices into which tail is divided
  cl0 <- cl # Initial confidence level
  term <- VaR
  
  delta.cl <- (1 - cl)/number.slices # Increment to confidence level as each slice is taken
  for (i in 1:(number.slices - 1)) {
    cl <- cl0 + i * delta.cl # Revised cl
    term <- term + mu * matrix(1, dim(cl)[1], dim(cl)[2]) - sigma * log(-n * log(cl)) # NB Gumber term
  }
  v <- term / (number.slices - 1)
  
  # Plotting
  plot(cl0, v, xlab = "Confidence Level", ylab = "ES", type = "l")
  legend("topleft",c('Input parameters',
                     paste('Location parameter for daily L/P = ', mu),
                     paste('Scale parameter for daily L/P = ', sigma),
                     paste('Holding period = ', hp, ' days')), bty="n")
  title("Gumbel ES against confidence level")
  
} 