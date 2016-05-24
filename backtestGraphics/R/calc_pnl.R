#' Calculate cumulative pnl, mean pnl, dollar sharpe ratio
#' 
#' This function takes in a data frame and calculates the cumulative P&L, daily
#' P&L, annualized P&L, P&L volatility and dollar sharpe ratio as well as three
#' biggest drawdowns. There are also average daily return rate and the standard 
#' deviation of daily return rates.
#'  
#' @param x A data frame that contains data for individual commodities.
#' @param trade.freq The trading frequency of the portfolio, numeric
#' 
#' @return a list ith cumulative pnl, mean pnl, annualized volatility of pnl, 
#' dollar sharpe ratio and drawdown

calc_pnl <- function(x, trade.freq = 7){
  
  ## Make a list object to store all outputs
  
  f <- list()
  
  ## k is the multiple for annualized pnl 
  
  k <- ifelse(trade.freq == 1, sqrt(252), floor(sqrt(365 / trade.freq)))
  
  ## Calculate the cumulative pnl
  
  f$pnl.cum <- sum(x[["pnl"]], na.rm = TRUE)
  
  ## Calculate the mean pnl
  
  f$pnl.mean <- f$pnl.cum / nrow(x)
  
  ## Calculate the annualized pnl
  
  f$pnl.annualized <- f$pnl.mean * k^2
  
  ## Calculate the annualized volatility of pnl
  
  f$pnl.vol <- sd(x[["pnl"]], na.rm = TRUE) * k
  
  ## Calculate the dollar sharpe ratio 
  
  f$pnl.sharpe <- f$pnl.mean*k/sd(x[["pnl"]], na.rm = TRUE)
  
  ## Calculate the three biggest drawdowns
  
  f$pnl.drawdown <- drawdown(x)
  
  ## Calculate the mean return rate by averaging all daily return rates
  
  f$ret.mean <- mean(x[["ret"]], na.rm = TRUE)
  
  ## Calculate the standard deviation of return rates
  
  f$ret.sd <- sd(x[["ret"]], na.rm = TRUE)
  
  return(f)
}
