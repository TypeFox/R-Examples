#' Calculate Statistics
#' 
#' This helper function calculates the summary of the input data set. The 
#' function takes in the sliced data set and the calculated trading frequency,
#' and calculate statistics about the market value, returns and volatility for
#' later display.
#' 
#' @param x.list is the input data set. Such data set has been sliced by
#' slice_data function, and is labeled as "intermediate data set".
#' 
#' @return f A list that contains all calculated numbers. 

stat_calculation <- function(x.list){
  
  ## Initiate a list to contain all the summary statistics about the data set
  
  f <- list()
  
  ## Calculate the first and last days
  
  f$day1 <- min(as.Date(x.list$x[["date"]], origin = "1970-01-01"))
  f$day2 <- max(as.Date(x.list$x[["date"]], origin = "1970-01-01"))
  
  ## Calculate pnl, sharpe ratio and best performers. There are too many things
  ## to be calculated about P&L, so we put the calculation in an individual
  ## function.
  
  f$pnl <- calc_pnl(x.list$x, x.list$trade.freq)
  
  ## Calculate the digits for rounding. If the data set contains real dollar 
  ## values as data, f$round = 0 and the backtestGraphics function will round
  ## away all the decimal numbers. If the data set contains weights asdata, the 
  ## backtestGraphics function will round all market values and profits to their
  ## 6th decimal digit
  
  f$rounding <- ifelse(f$pnl$pnl.mean > 1, 0, 6)
  
  ## Find the best/worst month of performance. The function for this is long
  ## enough that the best/worst month calculation deserves an individual
  ## function
  
  f$performance <- best_worst_month(x.list$x)
  
  f$gmv.mean <- mean(x.list$x[["gmv"]], na.rm = TRUE)
  
  ## k is the multiple for annualized pnl 
  
  k <- ifelse(x.list$trade.freq == 1, sqrt(252), floor(sqrt(360 / x.list$trade.freq)))
  
  ## Calculate annualized return rate
  
  f$meanret <- f$pnl$ret.mean * 100
  f$annualizedret <- f$meanret * k^2
  
  ## Calculate annualized volatility of return 
  
  f$volret <- f$pnl$ret.sd * 100 * k
  
  ## Calculate the sharpe ratio
  
  f$sharpe <- f$annualizedret / f$volret
  
  return(f)
}