#' Bivariate Gaussian Copule VaR
#' 
#' Derives VaR using bivariate Gaussian copula with specified inputs 
#' for normal marginals.
#' 
#' @param mu1 Mean of Profit/Loss on first position 
#' @param mu2 Mean of Profit/Loss on second position
#' @param sigma1 Standard Deviation of Profit/Loss on first position
#' @param sigma2 Standard Deviation of Profit/Loss on second position
#' @param rho Correlation between Profit/Loss on two positions
#' @param number.steps.in.copula Number of steps used in the copula approximation
#' ( approximation being needed because Gaussian copula lacks a closed form solution)
#' @param cl VaR confidece level
#' @return Copula based VaR
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Dowd, K. and Fackler, P. Estimating VaR with copulas. Financial Engineering
#' News, 2004.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # VaR using bivariate Gaussian for X and Y with given parameters:
#'    GaussianCopulaVaR(2.3, 4.1, 1.2, 1.5, .6, 10, .95)
#'
#' @export
GaussianCopulaVaR <- function(mu1, mu2, sigma1, sigma2, rho, 
                              number.steps.in.copula, cl){
  
  p <- 1 - cl # p is tail probability or cdf
  
  # For comparision, compute analytical normal VaR
  portfolio.mu <- mu1 + mu2
  portfolio.variance <- sigma1^2 + 2 * rho * sigma1 * sigma2 + sigma2^2
  portfolio.sigma <- sqrt(portfolio.variance)
  analytical.variance.covariance.VaR <- - portfolio.mu - portfolio.sigma * qnorm(p, 0, 1)
  # Specify bounds arbitrarily (NB: Would need to change manually if these were
  # inappropriate)
  L <- -portfolio.mu - 5 * portfolio.sigma
  fL <- CdfOfSumUsingGaussianCopula(L, mu1, mu2, sigma1, sigma2, 
                                    rho, number.steps.in.copula) - p
  sign.fL <- sign(fL)
  U <- -portfolio.mu + 5 *  portfolio.sigma
  fU <- CdfOfSumUsingGaussianCopula(U, mu1, mu2, sigma1, sigma2, 
                                    rho, number.steps.in.copula) - p
  sign.fU <- sign(fU)
  
  if (sign.fL == sign.fU){
    stop("Assumed bounds do not include answer")
  }
  
  # Bisection Algorithm
  tol <- 0.001 # Tolerance level (NM: change manually if desired)
  while (U - L > tol){
    x <- (L + U) / 2 # Bisection carried out in terms of P/L quantiles or minus VaR
    cum.prob <- CdfOfSumUsingGaussianCopula(x, mu1, mu2, sigma1, sigma2, rho, 
                                            number.steps.in.copula)
    fx <- cum.prob - p
    if (sign(fx) == sign(fL)){
      L <- x
      fL <- fx
    } else {
      U <- x
      fU <- fx
    }
  }
  
  y <- -x # VaR is negative of terminal x-value or P/L quantile
  
  cat("Analytical Variance Covariance Var: ", analytical.variance.covariance.VaR, "\n")
  cat("Var using bivariate gaussian copula: ", y, "\n")
  cat("Error in Copula VaR Estimate: ", analytical.variance.covariance.VaR-y)
  return(y)
}
