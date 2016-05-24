#' Bootstrapped VaR Confidence Interval
#'
#' Estimates the 90% confidence interval for bootstrapped VaR, for confidence
#' level and holding period implied by data frequency.
#'
#' @param Ra Vector corresponding to profit and loss distribution
#' @param number.resamples Number of samples to be taken in bootstrap procedure
#' @param cl Number corresponding to Value at Risk confidence level
#' @return 90% Confidence interval for bootstrapped VaR
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # To be modified with appropriate data.
#'    # Estimates 90% confidence interval for bootstrapped Var for 95% 
#'    # confidence interval 
#'    Ra <- rnorm(1000)
#'    BootstrapVaRConfInterval(Ra, 500, 0.95)
#'
#' @import bootstrap
#'
#' @export
BootstrapVaRConfInterval <- function(Ra, number.resamples, cl){
  
  # Determine if there are three arguments
  if (nargs() < 3){
    stop("Too few arguments")
  }
  if (nargs() > 3){
    stop("Too many arguments")
  }
  
  profit.loss.data <- as.vector(Ra)
  
  # Preprocess data
  unsorted.loss.data <- -profit.loss.data # Derives L/P data from input P/L data
  losses.data <- sort(unsorted.loss.data) # Puts losses in ascending order
  n <- length(losses.data)
  
  # Check that inputs have correct dimensions
  if (is.vector(cl) & (length(cl) != 1) ) {
    stop("Confidence level must be a scalar")
  }
  if (length(number.resamples) != 1) {
    stop("Number of resamples must be a scalar")
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
  
  # VaR estimation
  VaR <- bootstrap(losses.data, number.resamples, HSVaR, cl)$thetastar
  y <- quantile(VaR, c(.05, .95))
  return(y)
  
}