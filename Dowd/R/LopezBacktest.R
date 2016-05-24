#' First (binomial) Lopez forecast evaluation backtest score measure
#'
#' Derives the first Lopez (i.e. binomial) forecast evaluation score
#' for a VaR risk measurement model.
#' 
#' @param Ra Vector of portfolio of profit loss distribution
#' @param Rb Vector of corresponding VaR forecasts
#' @param cl VaR confidence level
#' @return Something
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Lopez, J. A. Methods for Evaluating Value-at-Risk Estimates. Federal Reserve
#' Bank of New York Economic Policy Review, 1998, p. 121.
#' 
#' Lopez, J. A. Regulatory Evaluations of Value-at-Risk Models. Journal of Risk
#' 1999, 37-64.
#' 
#' @author Dinesh Acharya
#' @examples
#'    
#'    # Has to be modified with appropriate data:
#'    # LopezBacktest for given parameters
#'    a <- rnorm(1*100)
#'    b <- abs(rnorm(1*100))+2
#'    LopezBacktest(a, b, 0.95)
#'
#' @export
LopezBacktest <- function(Ra, Rb, cl){
  
  profit.loss <- as.vector(Ra)
  VaR <- as.vector(Rb)
  
  n <- length(profit.loss)
  p <- 1-cl
  excess.loss <- profit.loss-VaR # Derives excess loss
  excess.loss <- excess.loss[excess.loss>0] # Gets rid of non-positive excess
  x <- length(excess.loss) # Score for each positive excess loss observation
  y <- (2/n)*sum(x-n*p)^2 # Score measure
  return(y)
  
}