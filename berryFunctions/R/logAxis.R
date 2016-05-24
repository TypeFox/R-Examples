#' Label logarithmic axes
#' 
#' Shortcut to calling \code{\link{logVals}}, \code{\link{axis}} and \code{\link{abline}}
#' 
#' @return An invisible list with 
#'        \item{vals}{Values for lines and label positions} 
#'        \item{labs}{Formatted values for labels} 
#'        \item{all}{Values for lines}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014
#' @seealso \code{\link{logVals}}, \code{\link{log10}}
#' @keywords aplot dplot
#' @export
#' @examples
#' 
#' x <- 10^runif(200, -1, 2)
#' plot(x, yaxt="n", log="y", pch=16)
#' logAxis(2)
#' # overplot vertical lines:
#' logAxis(2, expr=points(x, pch=16), base=1, col.axis=4, font=2)
#' 
#' # plots where log="x" is not possible:
#' hist(log10(x), breaks=20, col.axis="grey", main="")
#' logAxis(side=3, expr=hist(log10(x), breaks=20, add=TRUE, col=3))
#' # or just use the new logHist function (Feb 2016):
#' logHist(x, breaks=20, col=3)
#' 
#' # automatic calculation of from, to and base:
#' plot(1:3, axes=FALSE)
#' logAxis(1:2) # side can be a vector - nice, huh?
#' plot(-1:4, axes=FALSE)
#' logAxis(1:2) # treshold for base 1 instead of c(1,2,5) at 4 exponents exceeded.
#' 
#' plot(1:3, axes=FALSE)
#' logAxis(1:2, allticks=TRUE, lcol=NA)
#' 
#' par(mar=c(3,3,1,4))
#' plot(8:15) ;  logAxis(4) # with exponents if they are above 5
#' plot(10^(1:4), ylim=10^c(4,1), type="o", log="y") # reverse axis:
#' plot(10^(1:5), log="y"); logAxis(4, exponent=3) # different treshold
#' plot(10^(1:5), log="y"); logAxis(4, expon=3, base=c(1,2,5), expobase1=TRUE)
#' plot(-8:5); logAxis(4, allbase=c(1,2,5)) # In case you want to mislead...
#' 
#' @param side Which \code{\link{axis}} are to be labeled? Can be a vector within \code{1:4}. DEFAULT: 1
#' @param log Is the axis logarithmic by plot(log="x")? internal DEFAULT: \code{\link{par}("xlog")} or "ylog". DEFAULT: NULL
#' @param lcol Color of gridlines drawn in the graph with \code{\link{abline}}, NA to suppress. DEFAULT: "grey"
#' @param lty,lwd Type of gridlines. DEFAULT: 1
#' @param allticks Place all intermediate ticklines at the axis (without labelling). DEFAULT: FALSE
#' @param allargs List of arguments passed to axis for allticks=TRUE. DEFAULT: NULL
#' @param expr Expression drawing over the ablines, like (points(x,y). Can be code within {braces}.
#' @param las LabelAxisStyle for the orientation of the labels. DEFAULT: 1
#' @param from Lower exponent OR vector with data, as in \code{\link{logVals}}. DEFAULT based on \code{\link{par}("usr")}
#' @param to High end exponent. DEFAULT: internally based on par("usr")
#' @param Range Override from and to as range.
#' @param base Bases to be used in \code{\link{logVals}}. DEFAULT: c(1,2,5) or 1, depending on from and to.
#' @param big.mark Symbol separating thousands, eg. space, comma, dot, etc. see "format" and "prettyNum". DEFAULT: "'"
#' @param decimal.mark Character separating comma values, see "format" and "prettyNum". DEFAULT: "."
#' @param scientific See \code{\link{format}}. DEFAULT: FALSE
#' @param exponent Starting at which exponent should \code{\link{logVals}} return an expression with exponents? DEFAULT: 5
#' @param expobase1 Should "n * " be appended before 10^exp if n=1? DEFAULT: FALSE
#' @param allbase base for \code{$all} (for horizontal lines). DEFAULT: 1:9
#' @param box Draw box at the end to overplot \code{\link{abline}s} crossing the box? DEFAULT: TRUE
#' @param \dots Further arguments passed to axis, like \code{lwd, col.ticks, hadj, lty}, ...
#' 
logAxis <- function(
side=1,
log=NULL,
lcol="grey",
lty=1, lwd=1,
allticks=FALSE,
allargs=NULL,
expr,
las=1,
from,
to,
Range,
base,
big.mark="'",
decimal.mark=".",
scientific=FALSE,
exponent=5,
expobase1=FALSE,
allbase=1:9,
box=TRUE,
...)
{
for(side_i in side)
{ # loop around each side
if(side_i==1 | side_i==3) # vertical lines, labels at x-axis:
  {
  # set from and to:
  from_i <- if(missing(from)) par("usr")[1] else from
  to_i   <- if(missing(to))   par("usr")[2] else to
  # determine base:
  base_i <- if(missing(base)) {if(abs(to_i-from_i)>4) 1 else c(1,2,5)} else base
  # get labels and positions:
  lv <- logVals(from=from_i, to=to_i, Range=Range, base=base_i, big.mark=big.mark,
            decimal.mark=decimal.mark, scientific=scientific, exponent=exponent,
            expobase1=expobase1, allbase=allbase)
  # draw lines
  if(is.null(log)) log <- par("xlog")
  if(log) abline(v=lv$all,     col=lcol, lty=lty, lwd=lwd)
  else abline(v=log10(lv$all), col=lcol, lty=lty, lwd=lwd)
  }
else # horizontal lines, labels at y-axis:
  {
  # set from and to:
  from_i <- if(missing(from)) par("usr")[3] else from
  to_i   <- if(missing(to))   par("usr")[4] else to
  # determine base:
  base_i <- if(missing(base)) {if(abs(to_i-from_i)>4) 1 else c(1,2,5)} else base
  # get labels and positions:
  lv <- logVals(from=from_i, to=to_i, Range=Range, base=base_i, big.mark=big.mark,
            decimal.mark=decimal.mark, scientific=scientific, exponent=exponent,
            expobase1=expobase1, allbase=allbase)
  # draw lines
  if(is.null(log)) log <- par("ylog")
  if(log) abline(h=lv$all,     col=lcol, lty=lty, lwd=lwd)
  else abline(h=log10(lv$all), col=lcol, lty=lty, lwd=lwd)
  }
# axis labels:
if(log) axis(side=side_i, at=lv$vals,        labels=lv$labs, las=las, ...)
else    axis(side=side_i, at=log10(lv$vals), labels=lv$labs, las=las, ...)
if(allticks)
  {
  lv$add <- lv$all[!lv$all %in% lv$vals]
  if(!log) lv$add <- log10(lv$add)
  do.call(axis, args=owa(list(side=side_i, at=lv$add, labels=FALSE, col="grey", ...), allargs))
  }
} # End of loop
# Box to cover up the lines plotted over the existing box:
if(box) graphics::box("plot")
# overplot ablines with expr:
if(!missing(expr)) expr
# output:
return(invisible(lv))
} # End of function
