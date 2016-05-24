#' Derives prob ( X + Y < quantile) using Gaussian copula
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
#' @param rho Correlation between P/Ls on two positions
#' @param number.steps.in.copula The number of steps used in the copula approximation
#' @return Probability of X + Y being less than quantile
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Dowd, K. and Fackler, P. Estimating VaR with copulas. Financial Engineering
#' News, 2004.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Prob ( X + Y < q ) using Gaussian Copula for X with mean 2.3 and std. .2
#'    # and Y with mean 4.5 and std. 1.5 with beta 1.2 at 0.9 quantile
#'    CdfOfSumUsingGaussianCopula(0.9, 2.3, 4.5, 1.2, 1.5, 0.6, 15)
#'
#' @export
CdfOfSumUsingGaussianCopula <- function(quantile, mu1, mu2, sigma1, sigma2, rho, number.steps.in.copula){
  
  # Define w variable
  w.min <- 0.01
  w.max <- 0.99 # Specify min and max of w
  dw <- 0.01 # Specify incremental change in w
  w <- seq(w.min, w.max, dw) # Define w-variable going from min to max in units of size dw
  
  # Obtain 'first' and 'second' sets of copula values corresponding to 
  # values of w (NB: The difference between the 'first' and 'second' goes 
  # into the calculation of the approximate copula differential)
  first.copula <- double(length(w)) # initialize variables
  second.copula <- double(length(w))
  approximate.copula.differential <- double(length(w))
  for (i in 1:length(w)) {
    first.copula[i] <- GaussianCopula(w[i], 
                                    pnorm(quantile - qnorm(w[i], mu1, sigma1), mu2, sigma2), 
                                    0, number.steps.in.copula)
  }
  for (i in 2:length(w)) {
    second.copula[i] <- GaussianCopula(w[i] - dw, 
                                    pnorm(quantile - qnorm(w[i], mu1, sigma1), mu2, sigma2), 
                                    0, number.steps.in.copula)
  }
  
  # Obtain approximate copula differentials from the above
  for (i in 2:length(w)){
    approximate.copula.differential[i] <- (first.copula[i]-second.copula[i])/dw
  }
  
  y<- sum(approximate.copula.differential) * dw
  
  return(y)
}

GaussianCopula <- function(u, v, rho, number.steps){
  # Derives value of Gaussian Copula
  # u is value of first marginal for random variable X
  # v is value of second marginal for random variable Y
  # rho is correlation coefficient between X and Y
  # number.steps is used for each of x and y
  # NB: Applies to bivariate Gaussian copula only
  
  # Define x and y
  # NB: These are hypothetical domains of the cdfs; without loss of generality
  # we can take these to be standard normal percentiles; this is because we are
  # only using these to determine double integrals with limits u and v
  x.min <- -3
  x.max <- qnorm(u, 0, 1)
  dx <- (x.max - x.min)/number.steps
  x <- seq(x.min, x.max, dx)
  y.min <- -3
  y.max <- qnorm(v, 0, 1)
  dy <- (y.max - y.min)/number.steps
  y <- seq(y.min, y.max, dy)
  
  # Obtain copula as approximate double sum, which appriximates the 'true'
  # double integral
  term <- matrix(0, length(x),length(y))
  for (i in 1:length(x)) {
    for (j in 1:length(y)) {
      term[i, j] <- exp(-(x[i]^2-2*rho*x[i]*y[j]+y[j]^2)/(2*(1-rho^2)))
    }
  }
  
  y <- sum(term) * dx * dy / (2 * pi * (1 - rho^2)^0.5)
    
  return(y)
}
