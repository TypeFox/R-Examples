#' Blanco-Ihle forecast evaluation backtest measure
#'
#' Derives the Blanco-Ihle forecast evaluation loss measure for a VaR
#' risk measurement model.
#' 
#' @param Ra Vector of a portfolio profit and loss
#' @param Rb Vector of corresponding VaR forecasts
#' @param Rc Vector of corresponding Expected Tailed Loss forecasts
#' @param cl VaR confidence interval
#' @return First Blanco-Ihle score measure.
#' 
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' Blanco, C. and Ihle, G. How Good is Your Var? Using Backtesting to Assess
#' System Performance. Financial Engineering News, 1999.
#' 
#' @author Dinesh Acharya
#' @examples
#'    
#'    # Blanco-Ihle Backtest For Independence for given confidence level.
#'    # The VaR and ES are randomly generated.
#'    a <- rnorm(1*100)
#'    b <- abs(rnorm(1*100))+2
#'    c <- abs(rnorm(1*100))+2
#'    BlancoIhleBacktest(a, b, c, 0.95)
#'
#' @export
BlancoIhleBacktest <- function(Ra, Rb, Rc, cl){
  
  profit.loss <- as.vector(Ra)
  VaR <- as.vector(Rb)
  ETL <- as.vector(Rc)
  
  n <- length(profit.loss)
  p <- 1-cl
  excess.loss <- -profit.loss[-profit.loss>VaR] # Derives excess loss
  m <- length(excess.loss)
  
  benchmark <- double(m)
  score <- double(m)
  for (i in 1:m){
	benchmark[i] <- (ETL[i]-VaR[i])/VaR[i]
	score[i] <- (excess.loss[i]-VaR[i])/VaR[i]-benchmark[i]
  }
  
  # First Blanco-Ihle score measure
  y <- (2/n)*sum(score)^2
  return(y)
}