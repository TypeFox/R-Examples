#' Compare beta distributions
#' 
#' Visually understand the effect of the beta distribution parameters
#' 
#' @return None. Used for plotting.
#' @note Tries to find suitable subplot for axis labels. This works only for increasing parameter values.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2015
#' @seealso \code{\link{betaPlot}}
#' @keywords hplot distribution
#' @export
#' @examples
#' 
#' betaPlotComp()
#' betaPlotComp(oma=c(2,2,2,2), ylim=lim0(5.5), textargs=list(y=NA))
#' betaPlotComp(shape1=c(3,10,34), shape2=c(7,9,24))
#' 
#' @param shape1 Vector of alpha values as in \code{\link{dbeta}}. DEFAULT: c(0.5, 1:4, 10,20)
#' @param shape2 Beta values to be compared. DEFAULT: shape1
#' @param cumulative Should the cumulative density distribution line be added? DEFAULT: FALSE
#' @param cex Character EXpansion size. DEFAULT: 0.8
#' @param las Label Axis Style passed to \code{\link{axis}}. DEFAULT: 1
#' @param main Main as in \code{\link{plot}}. DEFAULT: ""
#' @param ylim LIMit for the Y axis. DEFAULT: lim0(4)
#' @param mar MARgins for plot passed to \code{\link{par}}. DEFAULT: rep(0,4)
#' @param oma Outer MArgins for plot passed to \code{\link{par}}. DEFAULT: c(2,2,4.5,2)
#' @param mgp MarGin Placement. DEFAULT: c(3,0.7,0)
#' @param keeppar Should margin parameters be kept instead of being restored to previous value? DEFAULT: FALSE
#' @param textargs List of arguments passed to \code{\link{textField}}. DEFAULT: NULL
#' @param \dots Further arguments passed to \code{\link{betaPlot}} like lines, fill, etc.
#' 
betaPlotComp <- function(
  shape1=c(0.5, 1:4, 10,20),
  shape2=shape1, 
  cumulative=FALSE,
  cex=0.8, 
  las=1, 
  main="", 
  ylim=lim0(4), 
  mar=rep(0,4), 
  oma=c(2,2,4.5,2),  
  mgp=c(3,0.7,0), 
  keeppar=FALSE,
  textargs=NULL,  
  ... 
  )
{
a <- shape1
b <- shape2
op <- par(mfrow=c(length(b), length(a)), yaxt="n", xaxt="n", cex=cex, mar=mar, oma=oma, mgp=mgp)
on.exit(if(!keeppar) par(op))
for(y in b)
  for(x in a)
  {
  betaPlot(x,y, main=main, ylim=ylim, mar=mar, cumulative=cumulative, ...)
  #mtext(text=paste("a=",x,", b=",y, sep=""), side=3, line=-1.5, cex=0.7)
  textdef <- list(x=0.5, y=0.9*ylim[2], labels=paste("a=",x,", b=",y, sep=""), 
                  mar=0.1, fill=addAlpha("white", 0.7), cex=1)
  do.call(textField, owa(textdef, textargs))
  if(y==min(b)) mtext(paste("a =",x), side=3, line=0.5, cex=cex)
  if(x==min(a)) mtext(paste("b =",y), side=2, line=0.5, cex=cex)
  if(y==max(b) & x==quantile(a, 0.5, type=1)) axis(1, at=c(0,0.5,1), xaxt="s", las=las)
  if(x==max(a) & y==quantile(b, 0.5, type=1)) axis(4, yaxt="s", las=las)
  }
mtext("beta density distribution", line=2.5, outer=TRUE, cex=1.2, font=2)
#box("outer")
}
