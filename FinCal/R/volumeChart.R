#' Technical analysts - Volume charts: show each period's volume as a vertical line
#' 
#' @param ohlc output from get.ohlc.yahoo or get.ohlc.google
#' @param main an overall title for the plot
#' @param ... Arguments to be passed to ggplot
#' @seealso \code{\link{get.ohlc.yahoo}}
#' @seealso \code{\link{get.ohlc.google}}
#' @importFrom ggplot2 ggplot geom_bar theme labs aes_string element_text
#' @export
#' @examples
#' # google <- get.ohlc.yahoo("GOOG");
#' # volumeChart(google)
#' # apple <- get.ohlc.google("AAPL");
#' # volumeChart(apple)
volumeChart <- function(ohlc,main="",...){
  options(warn=-1)
  ohlc$date <- as.Date(ohlc$date,"%Y-%m-%d")
  g=ggplot(ohlc, aes_string(x='date',y='volume')) + geom_bar(stat="identity") + labs( title =main ) + labs(x="") + 
    theme(axis.text.x=element_text(angle=90))
  return(g)
}
