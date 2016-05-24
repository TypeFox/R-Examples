#' Plot of Emperical and 2 Generalised Pareto mean excess functions 
#'
#' Plots of emperical mean excess function and two generalized pareto mean excess functions which differ in their tail-index value.
#'
#' @param Ra Vector of daily Profit/Loss data 
#' @param mu Location parameter
#' @param beta Scale parameter
#' @param zeta1 Assumed tail index for first mean excess function
#' @param zeta2 Assumed tail index for second mean excess function
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes ES assuming generalised Pareto for following parameters
#'    Ra <- 5 * rnorm(100)
#'    mu <- 1
#'    beta <- 1.2
#'    zeta1 <- 1.6
#'    zeta2 <- 2.2
#'    GParetoMultipleMEFPlot(Ra, mu, beta, zeta1, zeta2)
#'
#' @export
GParetoMultipleMEFPlot <- function(Ra, mu, beta, zeta1, zeta2) {
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
  gpmef1 <-  (1 + zeta1 * (u - mu) / beta)/(1 - zeta1);
  gpmef2 <-  (1 + zeta2 * (u - mu) / beta)/(1 - zeta2);
  # Plot
  # Limits of axis
  xlims <- c(min(u),max(u))
  ylims <- c(min(mef, gpmef1, gpmef2), max(mef, gpmef1, gpmef2))
  plot(u , mef, xlims, ylims, type = "l", xlab = "Threshold (u)", 
       col = 5, ylab = "e(u)")
  par(new = TRUE)
  plot(u , gpmef1, xlims, ylims, type = "l", xlab = "Threshold (u)", 
       col = 4, ylab = "e(u)")
  par(new = TRUE)
  plot(u , gpmef2, xlims, ylims, type = "l", xlab = "Threshold (u)", 
       col = 3, ylab = "e(u)")
  title("Emperical and Two Generalised Pareto MEFs")
  legend("topright", legend = c("Emperical MEF", "Generalized Pareto MEF1", "Generalized Pareto MEF1"), text.col = c(5,4,3))
  
}