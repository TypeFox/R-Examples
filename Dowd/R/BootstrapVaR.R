#' Bootstrapped VaR for specified confidence level
#'
#' Estimates the bootstrapped VaR for confidence level and holding period
#' implied by data frequency.
#'
#' @param Ra Vector corresponding to profit and loss distribution
#' @param number.resamples Number of samples to be taken in bootstrap procedure
#' @param cl Number corresponding to Value at Risk confidence level
#' @return Bootstrapped VaR
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates bootstrapped VaR for given parameters
#'    a <- rnorm(100) # generate a random profit/loss vector
#'    BootstrapES(a, 50, 0.95)
#'
#' @import bootstrap
#'
#' @export
BootstrapVaR <- function(Ra, number.resamples, cl){
  
  if (nargs() < 3){
    stop("Too few arguments")
  }
  if (nargs() > 3){
    stop("Too many arguments")
  }
  
  profit.loss.data <- as.vector(Ra)
  # Preprocess data
  unsorted.loss.data <- -profit.loss.data
  losses.data <- sort(unsorted.loss.data)
  n <- length(losses.data)
  
  # Check that inputs have correct dimensions
  if (length(cl) != 1) {
    stop("Confidence level must be a scalar")
  }
  if (length(number.resamples) != 1){
    stop("Number of resamples must be a scalar");
  }
  
  # Check that inputs obey sign and value restrictions
  if (cl >= 1){
    stop("Confidence level must be less that 1")
  }
  if (cl <= 0){
    stop("Confidence level must be at least 0")
  }
  if (number.resamples <= 0){
    stop("Number of resamples must be at least 0")
  }
  
  # ES estimation
  VaR <- bootstrap(losses.data, number.resamples, HSVaR, cl)$thetastar
  y <- mean(VaR)
  return (y)
}