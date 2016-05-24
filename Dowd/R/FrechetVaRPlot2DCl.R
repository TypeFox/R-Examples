#' Plots Frechet Value at Risk against Cl
#'
#' Plots the VaR of a portfolio against confidence level assuming extreme losses
#' are Frechet distributed, for specified range of confidence level and a given 
#' holding period.
#'
#' Note that the long-right-hand tail is fitted to losses, not profits.
#' 
#' @param mu Location parameter for daily L/P
#' @param sigma Scale parameter for daily L/P
#' @param tail.index Tail index
#' @param n Block size from which maxima are drawn
#' @param cl Confidence level and should be a vector
#' @param hp Holding period and should be a scalar
#' 
#' @references Dowd, K. Measurh  ing Market Risk, Wiley, 2007.
#' 
#' Embrechts, P., Kluppelberg, C. and Mikosch, T., Modelling Extremal Events for
#' Insurance and Finance. Springer, Berlin, 1997, p. 324.
#' 
#' Reiss, R. D. and Thomas, M. Statistical Analysis of Extreme Values from 
#' Insurance, Finance, Hydrology and Other Fields, Birkhaueser, Basel, 1997, 
#' 15-18.
#' 
#' @author Dinesh Acharya
#' 
#' @examples
#' 
#'    # Plots VaR against vector of cl assuming Frechet Distribution for given parameters
#'    cl <- seq(0.9, .99, .01)
#'    FrechetVaRPlot2DCl(3.5, 2.3, 1.6, 10, cl, 30)
#'
#' @export
FrechetVaRPlot2DCl <- function(mu, sigma, tail.index, n, cl, hp){
  
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
  
  # Change cl to row vector
  cl <- t(as.matrix(cl))
  
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
  # Plotting
  plot(cl, VaR, type = "l", xlab = "Confidence level", ylab = "VaR", main = "Frechet VaR against confidence level")
  text(mean(cl),
       max(VaR) - .1*(max(VaR) - min(VaR)),
       'Input parameters')
  text(mean(cl),
       max(VaR)-.2*(max(VaR)-min(VaR)),
       paste('Location parameter for daily L/P = ', mu))
  text(mean(cl),
       max(VaR) - .3 * (max(VaR) - min(VaR)),
       paste('Scale parameter for daily L/P = ', sigma))
  text(mean(cl),
       max(VaR) - .4 * (max(VaR) - min(VaR)),
       paste('Tail index = ', tail.index))
  text(mean(cl),
       max(VaR) - .5 * (max(VaR) - min(VaR)),
       paste('Holding period = ', hp, ' days'))  
  
} 