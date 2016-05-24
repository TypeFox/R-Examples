#' Christoffersen Backtest for Independence
#'
#' Carries out the Christoffersen backtest of independence for a VaR risk
#' measurement model, for specified VaR confidence level.
#' 
#' @param Ra Vector of portfolio profit and loss observations
#' @param Rb Vector of corresponding VaR forecasts
#' @param cl Confidence interval for
#' @return Probability that given the data set, the null hypothesis 
#' (i.e. independence) is correct.
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Christoffersen, P. Evaluating Interval Forecasts. International Economic 
#' Review, 39(4), 1992, 841-862.
#' 
#' @author Dinesh Acharya
#' @author Dinesh Acharya
#' @examples
#'    
#'    # Has to be modified with appropriate data:
#'    # Christoffersen Backtest For Independence for given parameters
#'    a <- rnorm(1*100)
#'    b <- abs(rnorm(1*100))+2
#'    ChristoffersenBacktestForIndependence(a, b, 0.95)
#'
#' @export
ChristoffersenBacktestForIndependence <- function(Ra, Rb, cl){
  
  profit.loss <- as.vector(Ra)
  VaR <- as.vector(Ra)
  
  n <- length(profit.loss)
  p <- 1-cl
  excess.loss <- -profit.loss-VaR # Derives excess loss
  excess.loss <- excess.loss[excess.loss>0] # Gets rid of negative or zeros
  ##########################################
  # There are mistakes in original code and needs to be addressed.
  # VaR <- VaR[excess.loss>0]
  ##########################################
  t00 <- 0
  t01 <- 0
  t10 <- 0
  t11 <- 0
  for (i in 2:length(excess.loss)){
	if(excess.loss[i]<=0){
	  if(excess.loss[i-1]<=0){
		t00 <- t00+1
	  } else {
	    t10 <- t10+1
	  }
	} else {
	  if(excess.loss[i-1]<=0){
	    t01 <- t01+1
	  } else {
	    t11 <- t11+1
	  }
	}
  }
  
  # Recover pie terms
  pie0 <- t01/(t00+t01)
  pie1 <- t11/(t10+t11)
  
  # Likelihood ratio test statistic
  LR=-2*log((((1-p)^(t00+t10))*(p^(t01+t11)))+2*log((((1-pie0)^t00))*(pie0^t01)*((1-pie1)^(t10))*pie1^t11))
  
  # Probability that null hypothesis (independence is correct)
  y <- 1-pchisq(LR,1)
  return(y)
  
}