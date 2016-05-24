#' shift and scale a vector
#' 
#' rescale a numeric vector: map values linearly onto a given range
#' 
#' @return numeric vector, rescaled onto output range
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jan 2016
#' @seealso \code{scales::rescale}
#' @references \url{http://stackoverflow.com/a/18303620}
#' @keywords manip
#' @export
#' @examples
#' 
#' rescale(10:15, 135, 200)
#' rescale(10:15, 200, 135)
#' rescale(10:15, to=c(1,5))
#' 
#' values <- rbeta(1e3, shape1=4, shape2=35)
#' hist(rescale(values, 135, 200), breaks=25, col=3)
#' 
#' @param x Numerical vector of values to be mapped to a given range
#' @param from output minimum. DEFAULT: 0
#' @param to output maximum. DEFAULT: 1
#' 
rescale <- function(
x,
from=0,
to=1
)
{
if(length(from)!=1) warning("from has length ", length(from))
if(length(to)!=1)   warning(  "to has length ", length(to))
  (x-min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE)-min(x, na.rm=TRUE)) * (to - from) + from
}

