#' Plot of Emperical and Generalised Pareto mean excess functions 
#'
#' Plots of emperical mean excess function and Generalized mean excess function.
#'
#' @param Ra Vector of daily Profit/Loss data 
#' @param mu Location parameter
#' @param beta Scale parameter
#' @param zeta Assumed tail index
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes ES assuming generalised Pareto for following parameters
#'    Ra <- 5 * rnorm(100)
#'    mu <- 0
#'    beta <- 1.2
#'    zeta <- 1.6
#'    GParetoMEFPlot(Ra, mu, beta, zeta)
#'
#' @export
GParetoMEFPlot <- function(Ra, mu, beta, zeta) {
  x <- as.vector(Ra)
  x <- sort(x) 
  u <- x
  n <- length(u)
  mef <- double(n - 1)
  
  for (i in 1:(n - 1)) {
    x <- x[which(x > u[i])]
    mef[i] <- mean(x) - u[i]
  }
  
  u <- t(u)
  u <- u[u!=max(u)]
  gpmef <-  (1 + zeta * (u - mu) / beta)/(1 - zeta);
  # Plot
  # Limits of axis
  xlims <- c(min(u),max(u))
  ylims <- c(min(mef, gpmef), max(mef, gpmef))
  plot(u , mef, xlims, ylims, type = "l", xlab = "Threshold (u)", 
       col = 6, ylab = "e(u)")
  par(new = TRUE)
  plot(u , gpmef, xlims, ylims, type = "l", xlab = "Threshold (u)", 
       col = 3, ylab = "e(u)")
  title("Emperical and Generalised Pareto Mean Excess Functions")
  legend("topright", legend = c("Emperical MEF", "Generalized Pareto MEF"), text.col = c(6,3))
  
}