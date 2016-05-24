#' Download stock prices from Yahoo Finance (open, high, low, close, volume, adjusted)
#' 
#' @param symbol symbol of stock, e.g. AAPL, GOOG, SPX
#' @param start start date, e.g., 2013-07-31
#' @param end end date, e.g., 2013-08-06
#' @param freq time interval, e.g., d:daily, w:weekly, m:monthly
#' @seealso \code{\link{get.ohlcs.yahoo}}
#' @seealso \code{\link{get.ohlc.google}}
#' @importFrom utils read.csv
#' @export
#' @examples
#' # get.ohlc.yahoo(symbol="AAPL")
#' # get.ohlc.yahoo(symbol="AAPL",start="2013-08-01",freq="d")
#' # get.ohlc.yahoo(symbol="AAPL",start="2013-07-01",end="2013-08-01",freq="w")
get.ohlc.yahoo <- function(symbol,start="firstDay",end="today",freq="d"){
  if(start == "firstDay"){
      URL=paste("http://ichart.finance.yahoo.com/table.csv?s=", symbol, sep="")
      dat <- read.csv(URL)
      dat=dat[order(dat$Date),]
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
        URL=paste("http://ichart.finance.yahoo.com/table.csv?s=",symbol,"&a=",a,"&b=",b,"&c=",c,"&d=",d,"&e=",e,"&f=",f,"&g=",freq, sep="")
        dat <- read.csv(URL)
        dat=dat[order(dat$Date),]
  } 
colnames(dat) <- c("date", "open", "high", "low", "close", "volume", "adjusted")
return(dat)
}
