#' Normal density plot
#' 
#' Nice plot of normal density distribution
#' 
#' @details This function finds some nice defaults for very quickly plotting a normal distribution by just specifying mean and sd.
#' 
#' @return None. Used for plotting.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, July 2014
#' @seealso \code{\link{betaPlot}}, \code{\link{dnorm}},
#'       \url{http://cran.r-project.org/package=denstrip},
#'       \url{http://cran.r-project.org/web/views/Distributions.html}
#' @keywords hplot distribution
#' @export
#' @examples
#' 
#' normPlot()
#' normPlot(81.7, 11.45)
#' normPlot(180, 11, quant=FALSE, width=2)
#' 
#' @param mean average value as in \code{\link{dnorm}}. DEFAULT: 0
#' @param sd standard deviation. DEFAULT: 1
#' @param width distance (in sd) from plot ends to mean. DEFAULT: 3
#' @param lines Should vertical lines be plotted at mean +- n*sd? DEFAULT: TRUE
#' @param quant should quantile regions be drawn with \code{fill} colors? DEFAULT: TRUE
#' @param fill color(s) passed to \code{\link{polygon}}. DEFAULT: addAlpha("blue",c(2:6,7:2)/10)
#' @param cumulative Should cumulative density distribution be added? DEFAULT: TRUE
#' @param las arguments passed to \code{\link{plot}}. DEFAULT: 1
#' @param main main as in \code{\link{plot}}. DEFAULT: paste("Normal density with\\nmean =", mean, "and sd =", sd)
#' @param ylim limit for the y axis. DEFAULT: lim0(y)
#' @param ylab,xlab labels for the axes. DEFAULT: ""
#' @param type,lty,col arguments passed to \code{\link{lines}}. type="l" to add pdf line
#' @param mar margins for plot passed to \code{\link{par}}. DEFAULT: c(2,3,3,3)
#' @param keeppar should margin parameters be kept instead of being restored to previous value? DEFAULT: FALSE
#' @param \dots further arguments passed to \code{\link{plot}} like lwd, xaxs, cex.axis, etc.
#' 
normPlot <- function(
mean=0,
sd=1,
width=3,
lines=TRUE,
quant=TRUE,
fill=addAlpha("blue",c(2:6,7:2)/10),
cumulative=TRUE,
las=1,
main=paste("Normal density with\nmean =", signif(mean,2), "and sd =", signif(sd,2)),
ylim=lim0(dnorm(mean,mean,sd)),
ylab="",
xlab="",
type="n",
lty=1,
col=par("fg"),
mar=c(2,3,3,3),
keeppar=FALSE,
...
)
{
op <- par(mar=mar)
# set up plot
xp <- c(mean-width*sd, mean+width*sd)
plot(xp, 0:1, ylim=ylim, las=las, type="n", ylab=ylab, xlab=xlab, main=main, ...)
# vertical lines at   mean +- n*sd
if(lines) abline(v=mean+floor(-width:width)*sd, col=8)
# recalculate coordinates with extrema, so that polygon goes to zero
x <- seq(par("usr")[1], par("usr")[2], length=400)
x <- c(mean-20*sd, x, mean+20*sd)
y <- dnorm(x, mean, sd)
# and plot that
if(!quant)
{
#fill <- fill[ceiling(length(fill)/2)] # get the middle color value
if(length(fill)>1) fill <- fill[4]
polygon(c(x,x[1]), c(0,y[-1],0), col=fill, border=NA)
}
else
{
qvals <- seq(0.05,0.95,by=0.1) # quantile borders
fill <- rep(fill, length=11) # fill recycling
xvals <- c(par("usr")[1], qnorm(qvals, mean, sd),par("usr")[2])
for(i in 2:length(xvals))
  {
  xv <- seq(xvals[i-1], xvals[i], length=30)
  yv <- dnorm(xv, mean, sd)
  polygon(c(xv,tail(xv,1),xv[1]), c(yv,0,0), col=fill[i-1], border=NA)
  }
yvals <- seq(ylim[2]/3, ylim[2]/20, len=5)
text(xvals[7:11], yvals, c(10,30,50,70,90), adj=c(1,-0.1))
segments(x0=xvals[6:2],x1=xvals[7:11],y0=yvals)
}
# PDF (probability density function) line:
lines(x, y, type=type, col=col, lty=lty)
# CDF (cumulated density function) line:
if(cumulative)
  {
  lines(x, pnorm(x, mean, sd)*par("usr")[4], type="l", col=2)
  axis(4, at=0:4/4*par("usr")[4], labels=0:4/4, col.axis=2, las=1, col=2)
  }
box()
if(!keeppar) par(op)
}
