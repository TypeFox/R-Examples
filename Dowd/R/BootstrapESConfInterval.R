#' Bootstrapped ES Confidence Interval
#'
#' Estimates the 90% confidence interval for bootstrapped ES, for confidence
#' level and holding period implied by data frequency.
#'
#' @param Ra Vector corresponding to profit and loss distribution
#' @param number.resamples Number of samples to be taken in bootstrap procedure
#' @param cl Number corresponding to Expected Shortfall confidence level
#' @return 90% Confidence interval for bootstrapped ES
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # To be modified with appropriate data.
#'    # Estimates 90% confidence interval for bootstrapped ES for 95% 
#'    # confidence interval 
#'    Ra <- rnorm(1000)
#'    BootstrapESConfInterval(Ra, 50, 0.95)
#'
#' @import bootstrap
#'
#' @export
BootstrapESConfInterval <- function(Ra, number.resamples, cl){
  
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
 
  # ES estimation
  es <- bootstrap(losses.data, number.resamples, HSES, cl)$thetastar
  y <- quantile(es, c(.05, .95))
  return(y)
  
}
