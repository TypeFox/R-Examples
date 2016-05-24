#' Find Trading Frequency
#' 
#' This helper function looks for the trading frequency of the input data set.
#' It can identify daily, weekly, monthly or yearly trading frequency in the 
#' data set. 
#' 
#' @param x The input data set whose trading period we are interested in
#' 
#' @return trade.freq The number that indicates the trading frequency. 

trade_freq <- function(x){
  
  date.unique <- sort(unique(x$date))
  
  ## Create two lists of dates by cuting the first and the last trading dates.
  ## Then subtract one list from the other to calculate the difference between
  ## any two neighboring trading dates.
  
  date1       <- date.unique[-length(date.unique)]
  date2       <- date.unique[-1]
  date.diff   <- date2 - date1
  
  ## Look at the differences between neighboring trading dates and decide
  ## trading periods. If trading period does not look like daily, weekly, 
  ## monthly nor yearly, the function returns an error
  
  if(length(which(date.diff <= 3)) > 0.3 * length(date.diff)){
    trade.freq <- 1
  }
  else if(length(which(date.diff >= 5 & date.diff <= 10)) > 0.3 * length(date.diff)){
    trade.freq <- 7
  }
  else if(length(which(date.diff >= 28 & date.diff <= 31)) > 0.3 * length(date.diff)){
    trade.freq <- 30
  }
  else if(length(which(date.diff >= 80 & date.diff <= 100)) > 0.3 * length(date.diff)){
    trade.freq <- 90
  }
  else{
    stop("It seems that the trading frequency of the input data set is not daily 
         or weekly or monthly or yearly. Please specify that manually.")
  }
  
  return(trade.freq)
}