#' @title Plots historical simulation ES against confidence level 
#'
#' @description Function plots the historical simulation ES of a 
#' portfolio against confidence level, for specified range of confidence level 
#' and holding period implied by data frequency.
#'
#' @param Ra Vector of daily P/L data
#' @param cl Vector of ES confidence levels
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots historical simulation ES against confidence level
#'    Ra <- rnorm(100)
#'    cl <- seq(.90, .99, .01)
#'    HSESPlot2DCl(Ra, cl)
#'
#' @export
HSESPlot2DCl <- function(Ra, cl){
  
  # Determine if there are three arguments, and ensure that arguments are read 
  # as intended
  if (nargs() < 2) {
    stop("Too few arguments.")
  }
  if (nargs() > 2) {
    stop("Too many arguments")
  }
  if (nargs() == 2) {
    profit.loss <- as.vector(Ra)
    n <- length(profit.loss)
  }
  
  # Check that inputs have correct dimensions
  cl <- as.matrix(cl)
  cl.rows <- dim(cl)[1]
  cl.columns <- dim(cl)[2]
  if (min(cl.rows, cl.columns) > 1) {
    stop("Confidence level must be a vector.")
  }
  cl <- as.vector(cl)
  
  # Check that inputs obey sign and value restrictions
  if (max(cl) >= 1) {
    stop("Confidence level must be less than 1.")
  }
  if (min(cl) <= 0) {
    stop("Confidence level must positive.")
  }
  
  # ES estimation
  ES <- double(length(cl))
  for (k in 1:length(cl)) {
    ES[k] <-  HSES(profit.loss, cl[k])
  }
  
  # Plot
  plot(cl, ES, type = "l", col = 5, xlab = "Confidence Level", ylab = "ES", main = "Historical ES against confidence level")
} 