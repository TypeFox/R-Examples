#' abline at locator point in graph
#'
#' Draw vertical and/or horizontal lines at positions in a graph located by clicking
#'
#' @details Not tested across platforms yet...
#' @return \code{\link{locator}} result
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Mar 2016
#' @seealso \code{\link{locator}}, \code{\link{abline}}
#' @keywords aplot iplot
#' @export
#' @examples
#'
#' plot(cumsum(rnorm(60)), type="l")
#' if(interactive()) locLine()

#' @param h Draw horizontal line at clicked location? DEFAULT: TRUE
#' @param v Draw vertical line at clicked location? DEFAULT: TRUE
#' @param n Number of points to be clicked. DEFAULT: 1
#' @param \dots Further arguments passed to \code{\link{abline}} like lty, lwd, col, etc

locLine <- function(
h=TRUE,
v=TRUE,
n=1,
...
)
{
coord <- locator(n=n)
if(h) abline(h=coord$y, ...)
if(v) abline(v=coord$x, ...)
coord
}
