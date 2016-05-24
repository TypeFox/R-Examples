#' nullplot
#' 
#' Make an plot with nothing in it
#' 
#' @param x1 lowest x-axis value
#' @param x2 largest x-axis value
#' @param y1 lowest y-axis value
#' @param y2 largest y-axis value
#' @param xlab x-axis title, defaults to no title
#' @param ylab y-axis title, defaults to no title
#' @param ... further arguments passed on to plot
#' @examples
#' nullplot()
nullplot <- function(x1=0,x2=1,y1=0,y2=1,xlab="",ylab="",...) {
  plot(0,0,xlim=c(x1,x2),ylim=c(y1,y2),type="n",xlab=xlab,ylab=ylab,...)
}
