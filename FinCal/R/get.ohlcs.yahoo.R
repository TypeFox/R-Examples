#' Batch download stock prices from Yahoo Finance (open, high, low, close, volume, adjusted)
#' 
#' @param symbols symbols of stock, e.g. AAPL, GOOG, SPX
#' @param start start date, e.g., 2013-07-31
#' @param end end date, e.g., 2013-08-06
#' @param freq time interval, e.g., d:daily, w:weekly, m:monthly
#' @seealso \code{\link{get.ohlc.yahoo}}
#' @seealso \code{\link{get.ohlcs.google}}
#' @importFrom utils read.csv
#' @export
#' @examples
#' # get.ohlcs.yahoo(symbols=c("AAPL","GOOG","SPY"),freq="d")
#' # get.ohlcs.yahoo(symbols=c("AAPL","GOOG","SPY"),start="2013-01-01",freq="m")
get.ohlcs.yahoo <- function(symbols,start="firstDay",end="today",freq="d"){
  n <- length(symbols)
  ohlc=list()
  if(start == "firstDay"){
      for(i in 1:n){
        URL=paste("http://ichart.finance.yahoo.com/table.csv?s=", symbols[i], sep="")
        dat <- read.csv(URL)
        dat=dat[order(dat$Date),]
        colnames(dat) <- c("date", "open", "high", "low", "close", "volume", "adjusted")
        ohlc[[symbols[i]]]=dat
      }
  }else{
    temp <- strsplit(start,"-")
    a=sprintf('%.2d',as.integer(temp[[1]][2])-1)
    b=temp[[1]][3]
    c=temp[[1]][1]
    if(end != "today"){
      temp <- strsplit(end,"-")
      d=sprintf('%.2d',as.integer(temp[[1]][2])-1)
      e=temp[[1]][3]
      f=temp[[1]][1]
    }else{
      end=as.character(Sys.Date())
      temp <- strsplit(end,"-")
      d=sprintf('%.2d',as.integer(temp[[1]][2])-1)
      e=temp[[1]][3]
      f=temp[[1]][1]
    }
    
    for(i in 1:n){
      URL=paste("http://ichart.finance.yahoo.com/table.csv?s=",symbols[i],"&a=",a,"&b=",b,"&c=",c,"&d=",d,"&e=",e,"&f=",f,"&g=",freq, sep="")
      dat <- read.csv(URL)
      dat=dat[order(dat$Date),]
      colnames(dat) <- c("date", "open", "high", "low", "close", "volume", "adjusted")
      ohlc[[symbols[i]]]=dat
    }
  } 
  return(ohlc)
}
