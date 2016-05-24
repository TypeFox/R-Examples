#' @title Plots historical simulation VaR and ES against confidence level 
#'
#' @description Function plots the historical simulation VaR and ES of a 
#' portfolio against confidence level, for specified range of confidence level 
#' and holding period implied by data frequency.
#'
#' @param Ra Vector of daily P/L data
#' @param cl Vectof of VaR confidence levels
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots historical simulation VaR and ES against confidence level
#'    Ra <- rnorm(100)
#'    cl <- seq(.90, .99, .01)
#'    HSVaRESPlot2DCl(Ra, cl)
#'
#' @export
HSVaRESPlot2DCl <- function(Ra, cl){
  
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
  
  # VaR and ES estimation
  VaR <- double(length(cl))
  ES <- double(length(cl))
  for (k in 1:length(cl)) {
    VaR[k] <- HSVaR(profit.loss, cl[k])
    ES[k] <-  HSES(profit.loss, cl[k])
  }
  
  # Plot
  plot(cl, VaR, type = "l", col = 3, xlab = "Confidence Level", ylab = "VaR/ES", main = "Historical VaR and ES against confidence level")
  par(new = TRUE)
  plot(cl, ES, type = "l", col = 5, xlab = "Confidence Level", ylab = "VaR/ES", main = "Historical VaR and ES against confidence level")
  legend("topleft", c("ES", "VaR"))
} 