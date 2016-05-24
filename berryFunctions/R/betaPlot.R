#' Beta density plot
#' 
#' Quick and nice plot of beta density distribution based on just alpha and beta
#' 
#' @details This function very quickly plots a beta distribution by just specifying alpha and beta.
#' 
#' @return None. Used for plotting.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, July 2014
#' @seealso \code{\link{betaPlotComp}}, \code{\link{normPlot}}, \code{\link{dbeta}}, 
#'    \url{http://cran.r-project.org/package=denstrip}, \url{http://cran.r-project.org/web/views/Distributions.html}
#' @keywords hplot distribution
#' @export
#' @examples
#' 
#' betaPlot()
#' betaPlot(2,1) 
#' betaPlot(0.5, 2)
#' 
#' # beta distribution is often used for proportions or probabilities
#' # overview of parameters
#' # alpha = number of successes + 1.  beta = number of failures + 1
#' betaPlotComp()
#' # a bigger: HDI (Highest Density Interval) further to the right (1)
#' # b bigger: HDI more to the left (0)
#' # both bigger: narrower HDI, stronger peak
#' 
#' @param shape1 Alpha value as in \code{\link{dbeta}}. DEFAULT: 1.5
#' @param shape2 Beta value. DEFAULT: 5
#' @param lines Quantiles at which vertical lines should be plotted. DEFAULT: NA
#' @param fill Color passed to \code{\link{polygon}}. DEFAULT: rgb(0,0.3,0.8, 0.4)
#' @param cumulative Should cumulative density distribution be added? DEFAULT: TRUE
#' @param mar Margins for plot passed to \code{\link{par}}. DEFAULT: c(2,3,3,3)
#' @param keeppar Should margin parameters be kept instead of being restored to previous value? DEFAULT: FALSE
#' @param las Label orientation, argument passed to \code{\link{plot}}. DEFAULT: 1
#' @param main main as in \code{\link{plot}}. DEFAULT: paste("Beta density with\\nalpha =", shape1, "and beta =", shape2)
#' @param ylim,xlim limit for the y and x axis. DEFAULT: lim0(y), 0:1
#' @param ylab,xlab labels for the axes. DEFAULT: ""
#' @param type,lty,col arguments passed to \code{\link{plot}} and \code{\link{lines}}.
#' @param \dots further arguments passed to \code{\link{plot}} like lwd, xaxs, cex.axis, etc.
#' 
betaPlot <- function(
  shape1=1.5,
  shape2=5, 
  lines=NA,
  fill=rgb(0,0.3,0.8, 0.4), 
  cumulative=TRUE,
  mar=c(2,3,3,3), 
  keeppar=FALSE,
  las=1, 
  main=paste("Beta density with\nalpha =", signif(shape1,3), "and beta =", signif(shape2,3)), 
  ylim=lim0(y),
  xlim=0:1,
  ylab="", 
  xlab="",
  type="l", 
  lty=1,
  col=par("fg"),
  ... 
  )
{
op <- par(mar=mar)
# create x and y coordinates
x <- seq(1e-3,1-1e-3, length=200)
y <- dbeta(x, shape1, shape2)
# plot
plot(x, y, las=las, type=type, col=col, lty=lty, xaxs="i", 
     ylab=ylab, xlab=xlab, main=main, ylim=ylim, xlim=xlim, ...)
# vertical lines
abline(v=qbeta(lines, shape1, shape2), col=8)
# Polygon
#browser()
polygon(c(0,x,1), c(0,y, 0), col=fill, border=NA)
lines(x, y, type=type, col=col, lty=lty)
if(cumulative)
  {
  lines(x, pbeta(x, shape1, shape2)*par("usr")[4], type="l", col=2)
  axis(4, at=0:4/4*par("usr")[4], labels=0:4/4, col.axis=2, las=1, col=2)
  }
box()
if(!keeppar) par(op)
}
