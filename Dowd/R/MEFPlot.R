#' Mean Excess Function Plot
#'
#' Plots mean-excess function values of the data set.
#'
#' @param Ra Vector data
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots mean-excess function values
#'    Ra <- rnorm(1000)
#'    MEFPlot(Ra)
#'
#' @export
MEFPlot <- function(Ra){
  data <- as.vector(Ra)
  if (!is.vector(data)) {
    stop("Input should be a vector data.")
  }
  data <- sort(data)
  u <- data
  n <- length(u)
  mef <- double(n-1)
  for (i in 1:(n-1)) {
    data <- data[which(data > u[i])]
    mef[i] <- mean(data) - u[i]
  }
  u <- u[!u == max(u)]
  plot(u, mef, type = "l", xlab = "Threshold", ylab = "e(u)", 
       main = "Empirical Mean Excess Function")
  
} 
