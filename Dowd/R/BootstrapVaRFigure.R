#' Plots figure of bootstrapped VaR
#'
#' Plots figure for the bootstrapped VaR, for confidence
#' level and holding period implied by data frequency.
#'
#' @param Ra Vector corresponding to profit and loss distribution
#' @param number.resamples Number of samples to be taken in bootstrap procedure
#' @param cl Number corresponding to Value at Risk confidence level
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # To be modified with appropriate data.
#'    # Estimates 90% confidence interval for bootstrapped VaR for 95% 
#'    # confidence interval 
#'    Ra <- rnorm(1000)
#'    BootstrapESFigure(Ra, 500, 0.95)
#'
#' @import bootstrap
#'
#' @export
BootstrapVaRFigure <- function(Ra, number.resamples, cl){
  
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
  VaR <- bootstrap(losses.data, number.resamples, HSVaR, cl)$thetastar
  mean.VaR <- mean(VaR)
  std.VaR <- sd(VaR)
  min.VaR <- min(VaR)
  max.VaR <- max(VaR)
  ninety.five.perc.conf.interval <- quantile(VaR, c(.05, .95))
  
  # Histogram
  cl.for.label <- 100*cl
  hist(VaR, 30, xlab="VaR", ylab="Frequency", main=paste("Bootstrapped Historical Simulation VaR at", cl, "% Confidence Level"))
  
}