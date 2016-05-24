#' Bivariate Gumbel Copule VaR
#' 
#' Derives VaR using bivariate Gumbel or logistic copula with specified inputs 
#' for normal marginals.
#' 
#' @param mu1 Mean of Profit/Loss on first position 
#' @param mu2 Mean of Profit/Loss on second position
#' @param sigma1 Standard Deviation of Profit/Loss on first position
#' @param sigma2 Standard Deviation of Profit/Loss on second position
#' @param beta Gumber copula parameter (greater than 1)
#' @param cl VaR onfidece level
#' @return Copula based VaR
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Dowd, K. and Fackler, P. Estimating VaR with copulas. Financial Engineering
#' News, 2004.
#'
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # VaR using bivariate Gumbel for X and Y with given parameters:
#'    GumbelCopulaVaR(1.1, 3.1, 1.2, 1.5, 1.1, .95)
#'
#' @export
GumbelCopulaVaR <- function(mu1, mu2, sigma1, sigma2, beta, cl){
  
  if (beta <= 1) {
    stop("Beta must be bigger than 1")
  }
  
  p <- 1 - cl # p is tail probability or cdf
  
  # Compute portfolio mean and sigma (NB: These are used here to help compute
  # initial bounds automatically)
  portfolio.mu <- mu1 + mu2
  portfolio.variance <- sigma1^2+sigma2^2
  portfolio.sigma <- sqrt(portfolio.variance)
  
  # Specify bounds arbitrarily (NB: Would need to change manually if these were
  # inappropriate)
  L <- -portfolio.mu - 5 * portfolio.sigma
  fL <- CdfOfSumUsingGumbelCopula(L, mu1, mu2, sigma1, sigma2, beta) - p
  sign.fL <- sign(fL)
  U <- -portfolio.mu + 5 *  portfolio.sigma
  fU <- CdfOfSumUsingGumbelCopula(U, mu1, mu2, sigma1, sigma2, beta) - p
  sign.fU <- sign(fU)
  if (sign.fL == sign.fU){
    stop("Assumed bounds do not include answer")
  }
  
  # Bisection Algorithm
  # The tolerance level in Dowd's original code was 0.0001. It was changed for test case
  # to take less time.
  tol <- 0.001 # Tolerance level (NM: change manually if desired)
  while (U - L > tol){
    x <- (L + U) / 2 # Bisection carried out in terms of P/L quantiles or minus VaR
    cum.prob <- CdfOfSumUsingGumbelCopula(x, mu1, mu2, sigma1, sigma2, beta)
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
  return(y)
}