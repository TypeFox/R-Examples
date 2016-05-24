#' Histogram of logarithmic values
#' 
#' Draw histogram of values on a logarithmic scale with nice axis labels
#' 
#' @return none
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Feb 2016
#' @seealso \code{\link{logAxis}}, \code{\link{hist}}
#' @keywords aplot dplot
#' @export
#' @examples
#' 
#' dat <- rbeta(1e4, 2, 18)*100
#' hist(dat, col="tan", breaks=50)
#' logHist(dat, col="tan", breaks=50)
#' 
#' @param x Vector of numerical values
#' @param logargs A list of arguments passed to \code{\link{logAxis}}. DEFAULT: NULL
#' @param main Title of graph, internally from x. DEFAULT: xmain
#' @param xlab X axis label. DEFAULT: xname
#' @param \dots further arguments passed to \code{\link{hist}} like breaks, col, ..., but not xaxt or add.
#' 
logHist <- function(
x,
logargs=NULL,
main=xmain,
xlab=xname,
...)
{
xname <- deparse(substitute(x))
xmain <- paste0("Histogram of log10(",xname,")")
hist(x=log10(x), ..., main=main, xlab=xlab, xaxt="n")
do.call(logAxis, owa(list(), logargs))
hist(x=log10(x), ..., add=TRUE)
}

