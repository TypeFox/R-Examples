#' @title Figure of Historical SImulation VaR and ES and histogram of L/P
#'
#' @description Plots figure showing the historical simulation VaR and ES and histogram
#'  of L/P for specified confidence level and holding period implied by data 
#'  frequency.
#'
#' @param Ra Vector of profit loss data
#' @param cl VaR confidence level
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots figure showing VaR and histogram of P/L data
#'    Ra <- rnorm(100)
#'    HSESFigure(Ra, .95)
#'
#' @export
HSESFigure<- function(Ra, cl){
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
  if ( cl <= 0){
    stop("Confidence level must be positive.")
  }
  # VaR and ES estimation
  VaR <- HSVaR(profit.loss, cl) # HS VaR
  ES <- HSES(profit.loss, cl) # HS ES
  # Histogram
  n <- hist(-profit.loss, main = "Historical Simulation VaR", breaks = 50, col = 4,
            xlab = "Loss(+) / Profit(-)", ylab = "Frequency")

  # Insert line showing VaR
  v <- c(0, .625 * max(n$counts)) # Coordinates for VaR line
  u <- VaR * matrix(1, length(v),1) # Coordinates for VaR line
  lines(u, v, type = "l", col = "black") # VaR line

  # Insert line showing ES
  w <- c(0, .45 * max(n$counts)) # Coordinates for ES
  z <- ES * matrix(1, length(w), 1) # Coordinates for ES line
  lines(z, w, type = "l", col = "black") # ES line
  cl.for.label <- 100 * cl # Input to confidence level label
  

  # VaR line label
  text(VaR, .75*max(n$counts), paste('VaR at', cl.for.label, "%Cl"))
  text(VaR, .65*max(n$counts), c(round(VaR, digits = 2)))
  
  # ES line label
  text(ES, .55*max(n$counts), c('ES ='))
  text(ES, .45*max(n$counts), paste(round(ES, digits = 2)))
  
} 
