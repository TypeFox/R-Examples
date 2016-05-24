#' Technical analysts - Line charts: show prices for each period as a continuous line for multiple stocks
#' 
#' @param ohlcs output from get.ohlc.yahoo.mult or get.ohlc.google.mult
#' @param y y coordinates: close, open, high, low or adjusted (yahoo data only)
#' @param main an overall title for the plot
#' @param ... Arguments to be passed to ggplot
#' @seealso \code{\link{get.ohlcs.yahoo}}
#' @seealso \code{\link{get.ohlcs.google}}
#' @seealso \code{\link{lineChart}}
#' @importFrom ggplot2 ggplot geom_line theme scale_y_continuous labs aes_string element_text element_blank
#' @importFrom reshape2 melt
#' @export
#' @examples
#' # googapple <- get.ohlcs.yahoo(c("GOOG","AAPL"),start="2013-01-01"); 
#' # lineChartMult(googapple)
#' # googapple <- get.ohlcs.google(c("GOOG","AAPL"),start="2013-01-01"); 
#' # lineChartMult(googapple)
lineChartMult <- function(ohlcs,y="close",main="",...){
  options(warn=-1)
  sname=names(ohlcs)
  n=length(sname)
  Date=ohlcs[[sname[1]]]$date <- as.Date(ohlcs[[sname[1]]]$date,"%Y-%m-%d")
  z=data.frame(Date)
  colname=c("Date")
  for(i in 1:n){
    colname=append(colname,sname[i])
    if(y == "close"){
      z[,i+1] = ohlcs[[sname[i]]]$close
    }else if(y == "open"){
      z[,i+1] = ohlcs[[sname[i]]]$open
    }else if(y == "high"){
      z[,i+1] = ohlcs[[sname[i]]]$high
    }else if(y == "low"){
      z[,i+1] = ohlcs[[sname[i]]]$low
    }else if(y == "adjusted"){
      z[,i+1] = ohlcs[[sname[i]]]$adjusted
    }
  }
  colnames(z) <- colname
  pdf <- melt(z, id="Date")
  ggplot(data=pdf, aes_string(x='Date', y='value', colour='variable'),...) + geom_line() +
    labs( title =main ) + labs(x="") + scale_y_continuous(name=y) +
    theme(axis.text.x=element_text(angle=90)) + theme(legend.title=element_blank())
}
