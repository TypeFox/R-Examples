#' Download stock prices from Google Finance (open, high, low, close, volume)
#' 
#' @param symbol symbol of stock, e.g. AAPL, GOOG, SPX
#' @param start start date, e.g., 2013-07-31
#' @param end end date, e.g., 2013-08-06
#' @seealso \code{\link{get.ohlc.yahoo}}
#' @seealso \code{\link{get.ohlcs.google}}
#' @importFrom RCurl getURL
#' @importFrom utils read.csv
#' @export
#' @examples
#' # get.ohlc.google(symbol="AAPL")
#' # get.ohlc.google(symbol="AAPL",start="2013-08-01")
#' # get.ohlc.google(symbol="AAPL",start="2013-07-01",end="2013-08-01")
get.ohlc.google <- function(symbol,start="2013-01-01",end="today"){
    temp <- strsplit(start,"-")
    a=month.abb[as.numeric(temp[[1]][2])]
    b=temp[[1]][3]
    c=temp[[1]][1]
    if(end != "today"){
      temp <- strsplit(end,"-")
      d=month.abb[as.numeric(temp[[1]][2])]
      e=temp[[1]][3]
      f=temp[[1]][1]
    }else{
      end=as.character(Sys.Date())
      temp <- strsplit(end,"-")
      d=month.abb[as.numeric(temp[[1]][2])]
      e=temp[[1]][3]
      f=temp[[1]][1]
    }
    URL=paste("https://www.google.com/finance/historical?q=",symbol,"&output=csv","&startdate=",a,"+",b,"+",c,"&enddate=",d,"+",e,"+",f, sep="")
    myCsv <- getURL(URL, ssl.verifypeer = FALSE)
    dat <- read.csv(textConnection(myCsv))
    colnames(dat) <- c("date", "open", "high", "low", "close", "volume")
    dates=as.character(dat$date)
    for(i in 1:length(dates)){
      tempdates=strsplit(dates[i],"-")
      tempdates[[1]][2]=match(tempdates[[1]][2],month.abb)
      tempd=paste(tempdates[[1]][1],tempdates[[1]][2],tempdates[[1]][3],sep="-")
      dates[i]=format(as.Date(tempd,"%d-%m-%y"),"%Y-%m-%d")
    }
    dat$date=dates
    dat=dat[order(dat$date),]
    return(dat)
}
