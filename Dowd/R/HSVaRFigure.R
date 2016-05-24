#' @title Figure of Historical SImulation VaR and histogram of L/P
#'
#' @description Plots figure showing the historical simulation VaR and histogram
#'  of L/P for specified confidence level and holding period implied by data 
#'  frequency.
#'
#' @param Ra Vector of profit loss data
#' @param cl ES confidence level
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots figure showing VaR and histogram of P/L data
#'    Ra <- rnorm(100)
#'    HSVaRFigure(Ra, .95)
#'
#' @export
HSVaRFigure<- function(Ra, cl){
  # Determine if there are two arguments and ensure that arguments are read as intended
  if (nargs() < 2) {
    stop("Too few arguments")
  }
  if (nargs() > 2){
    stop("Too many arguments")
  }
  if (nargs() == 2){
    profit.loss <- as.vector(Ra)
    n <- length(profit.loss)
  }
  
  # Check that inputs have correct dimensions
  if (length(cl) != 1) {
    stop("Confidence level must be a scalar")
  }
  
  if ( cl >= 1){
    stop("Confidence level must be less than 1")
  }
  
  # VaR estimation
  VaR <- HSVaR(profit.loss, cl) # HS VaR
  # Histogram
  n <- hist(profit.loss, main = "Historical Simulation VaR", col = 4,
            xlab = "Loss(+) / Profit(-)", ylab = "Frequency")
  v <- c(0, .625 * max(n$counts)) # Coordinates for VaR line
  u <- VaR * matrix(1, length(v),1) # Coordinates for VaR line
  lines(u, v, type = "l", col="blue")
  cl.for.label <- 100 * cl
  legend("topleft",c(paste('VaR at', cl.for.label, '% CL'), 
                      paste('=', VaR)), bty="n", cex = 0.7)
  
} 