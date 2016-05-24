#' Corn-Fisher ES
#'
#' Function estimates the ES for near-normal P/L using the Cornish-Fisher
#' adjustment for non-normality, for specified confidence level.
#'
#' @param mu Mean of P/L distribution
#' @param sigma Variance of of P/L distribution
#' @param skew Skew of P/L distribution
#' @param kurtosis Kurtosis of P/L distribution
#' @param cl ES confidence level
#' @return Expected Shortfall
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Zangri, P. A VaR methodology for portfolios that include options. 
#' RiskMetrics Monitor, First quarter, 1996, p. 4-12.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates Cornish-Fisher ES for given parameters
#'    CornishFisherES(3.2, 5.6, 2, 3, .9)
#'
#' @export
CornishFisherES <- function(mu, sigma, skew, kurtosis, cl) {
  # ES estimation
  confidence.level <- cl
  delta.cl <- (1 - confidence.level)/100
  cl <- double(99)
  z <- double(99)
  VaR <- double(99)
  adjustment <- double(99)
  cl[1] <- confidence.level
  for (i in 1:98){
    cl [i] <- cl[1] + i * delta.cl
    z[i] <- qnorm(1 - cl[i], 0, 1)
    adjustment[i] <- (1/6) * (z[i]^2 - 1) * skew + 
      (1/24) * (z[i]^3 - 3 * z[i]) * (kurtosis - 3) -
      (1/36) * (2 * z[i]^3 - 5 * z[i]) * skew^2
    VaR[i] <- - sigma * (z[i] + adjustment[i]) - mu
  }
  y <- mean(VaR)
  return(y)
  
}