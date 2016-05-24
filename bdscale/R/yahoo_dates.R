#' Get past trading days using close prices of supplied ticker
#' 
#' @param ticker The ticker you want to use, defaults to S&P 500: \code{^GSPC}
#' @return returns a vector of \code{Date}s
#' @importFrom utils URLencode read.csv
#' 
yahoo <- function(ticker='^GSPC') {
  encoded <- URLencode(ticker)
  u <- sprintf('http://real-chart.finance.yahoo.com/table.csv?s=%s&g=d&a=0&b=3&c=1950&ignore=.csv', encoded)
  csv <- read.csv(url(u))
  
  sort(as.Date(csv$Date))
}