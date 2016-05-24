#' Derives prob ( X + Y < quantile) using Gumbel copula
#' 
#' If X and Y are position P/Ls, then the VaR is equal to minus quantile. In
#' such cases, we insert the negative of the VaR as the quantile, and the
#' function gives us the value of 1 minus VaR confidence level. In other
#' words, if X and Y are position P/Ls, the quantile is the negative of the
#' VaR, and the output is 1 minus the VaR confidence level.
#' 
#' @param quantile Portfolio quantile (or negative of Var, if X, Y are position P/Ls)
#' @param mu1 Mean of Profit/Loss on first position 
#' @param mu2 Mean of Profit/Loss on second position
#' @param sigma1 Standard Deviation of Profit/Loss on first position
#' @param sigma2 Standard Deviation of Profit/Loss on second position
#' @param beta Gumber copula parameter (greater than 1)
#' @return Probability of X + Y being less than quantile
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Dowd, K. and Fackler, P. Estimating VaR with copulas. Financial Engineering
#' News, 2004.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Prob ( X + Y < q ) using Gumbel Copula for X with mean 2.3 and std. .2
#'    # and Y with mean 4.5 and std. 1.5 with beta 1.2 at 0.9 quantile
#'    CdfOfSumUsingGumbelCopula(0.9, 2.3, 4.5, 1.2, 1.5, 1.2)
#'
#' @export
CdfOfSumUsingGumbelCopula <- function(quantile, mu1, mu2, sigma1, sigma2, beta){
  
  if (beta <= 1) {
    stop("Beta must be bigger than 1")
  }
  
  # Define w variable
  w.min <- 0.001
  w.max <- 0.999 # Specify min and max of w
  dw <- 0.001 # Specify incremental change in w. Note that this is different than in Dowd's
  # original code (dw=0.001). It was necessary to reduce run time for examples.
  w <- seq(w.min, w.max, dw) # Define w-variable going from min to max in units of size dw
  
  # Obtain 'first' and 'second' sets of copula values corresponding to 
  # values of w (NB: The difference between the 'first' and 'second' goes 
  # into the calculation of the approximate copula differential)
  first.copula <- double(length(w)) # initialize variables
  second.copula <- double(length(w))
  approximate.copula.differential <- double(length(w))
  for (i in 1:length(w)) {
    first.copula[i] <- GumbelCopula(w[i], 
                                    pnorm(quantile - qnorm(w[i], mu1, sigma1), mu2, sigma2), 
                                    beta)
  }
  for (i in 2:length(w)) {
    second.copula[i] <- GumbelCopula(w[i] - dw, 
                                    pnorm(quantile - qnorm(w[i], mu1, sigma1), mu2, sigma2), 
                                    beta)
  }
  
  # Obtain approximate copula differentials from the above
  for (i in 2:length(w)){
    approximate.copula.differential[i] <- (first.copula[i]-second.copula[i])/dw
  }
  
  y<- sum(approximate.copula.differential) * dw
  
  return(y)
}

GumbelCopula <- function(u, v, beta){
  # Derives value of Gumbel Copula
  # u is value of first marginal for random variable X
  # v is value of second marginal for random variable Y
  # beta parameter for Gumbel Copula
  y <- exp(-((-log(u))^beta + (-log(v))^beta))^(1/beta)
  return(y)
}