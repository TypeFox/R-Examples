#' Corn-Fisher VaR
#'
#' Function estimates the VaR for near-normal P/L using the Cornish-Fisher
#' adjustment for non-normality, for specified confidence level.
#'
#' @param mu Mean of P/L distribution
#' @param sigma Variance of of P/L distribution
#' @param skew Skew of P/L distribution
#' @param kurtosis Kurtosis of P/L distribution
#' @param cl VaR confidence level
#' @return Value at Risk
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Zangri, P. A VaR methodology for portfolios that include options. 
#' RiskMetrics Monitor, First quarter, 1996, p. 4-12.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates Cornish-Fisher VaR for given parameters
#'    CornishFisherVaR(3.2, 5.6, 2, 3, .9)
#'
#' @export
CornishFisherVaR <- function(mu, sigma, skew, kurtosis, cl) {
  # Normal variate
  z <- qnorm(1 - cl, 0, 1)
  # Deriving adjustment factor
  adjustment <- (1/6) * (z^2 - 1) * skew + 
    (1/24) * (z^3 - 3 * z) * (kurtosis - 3) -
    (1/36) * (2 * z^3 - 5 * z) * skew^2
  # VaR estimation
  y <- - sigma * (z + adjustment) - mu
  return(y)
  
}