#' Plots figure of bootstrapped ES
#'
#' Plots figure for the bootstrapped ES, for confidence
#' level and holding period implied by data frequency.
#'
#' @param Ra Vector corresponding to profit and loss distribution
#' @param number.resamples Number of samples to be taken in bootstrap procedure
#' @param cl Number corresponding to Expected Shortfall confidence level
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
#'    BootstrapESFigure(Ra, 500, 0.95)
#'
#' @import bootstrap
#'
#' @export
BootstrapESFigure <- function(Ra, number.resamples, cl){
  
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
  
  # ES Estimation
  es <- bootstrap(losses.data, number.resamples, HSES, cl)$thetastar
  mean.es <- mean(es)
  std.es <- sd(es)
  min.es <- min(es)
  max.es <- max(es)
  ninety.five.perc.conf.interval <- quantile(es, c(.05, .95))
  
  # Histogram
  cl.for.label <- 100*cl
  hist(es, 30, xlab="ES", ylab="Frequency", main=paste("Bootstrapped Historical Simulation ES at", cl, "% Confidence Level"))
  
}