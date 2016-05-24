#' Carries out the binomial backtest for a VaR risk measurement model.
#' 
#' The basic idea behind binomial backtest (also called basic frequency 
#' test) is to test whether the observed frequency of losses that exceed VaR is
#' consistent with the frequency of tail losses predicted by the mode. Binomial
#' Backtest carries out the binomial backtest for a VaR risk measurement model 
#' for specified VaR confidence level and for a one-sided alternative 
#' hypothesis (H1).
#' 
#' @param x Number of failures
#' @param n Number of observations
#' @param cl Confidence level for VaR
#' @return Probability that the VaR model is correct
#' 
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' Kupiec, Paul. Techniques for verifying the accuracy of risk measurement
#' models, Journal of Derivatives, Winter 1995, p. 79. 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Probability that the VaR model is correct for 3 failures, 100 number
#'    # observations and  95% confidence level
#'    BinomialBacktest(55, 1000, 0.95)
#'
#' @export
BinomialBacktest <- function (x, n, cl){
  # Give warning if x>n or cl>100%
  if(x > n | cl >= 1 | cl <= 0){
    stop ("Incorrect parameter list. Make sure that x>n and 0<cl<1")
  }
  
  p <- 1 - cl # Probability of a failure each observation occurs
  
  if (x >= n*p){
    probability.model.is.correct <- 1-pbinom(x-1, n, p)
  } else {
    probability.model.is.correct <- pbinom(x, n, p)
  }
  return (probability.model.is.correct)
}
