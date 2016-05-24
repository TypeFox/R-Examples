#' Truncated pretty breakpoints
#' 
#' \code{\link{pretty}} with no values outside of x range
#' 
#' calculates \code{pretty(x)}, then removes the values that do not lie within
#' \code{\link{range}(x)}.\cr If force=TRUE, range(x) is reduced step by step
#' in a while loop until the condition is met. This is useful if you want
#' exactly 2 labels on an \code{\link{axis}}. In order not to get stuck, the
#' outer values are taken if there are more than n values within range(x).
#' 
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Aug 2014
#' @seealso \code{\link{pretty}} , \code{\link{logVals}}
#' @keywords dplot
#' @export
#' @examples
#' 
#' k <- c(135, 155, 120, 105, 140, 130, 190, 110)
#' range(k)
#' pretty(k)
#' pretty2(k)
#' 
#' pretty(c(0.2, 0.9), n=2)
#' pretty2(c(0.2, 0.9), n=2)
#' pretty2(c(0.2, 0.9), n=2, force=TRUE)
#' 
#' @param x object with numeric values
#' @param n desired number of values in \code{\link{pretty}}. DEFAULT: 5
#' @param force Must output lenght equal n exactly?  DEFAULT: FALSE
#' @param \dots all other arguments in \code{\link{pretty}}.
#' 
pretty2 <- function(
x,
n=5,
force=FALSE,
...)
{
p <- pretty(x, n=n, ...)
r <- range(x, finite=TRUE)
p <- p[p>=r[1] & p<=r[2]]
if(force)                   # Added March 2015
  while(length(p) != n)
  {
  x <- extendrange(x, f=-0.01)
  p <- pretty(x, n=n, ...)
  p <- p[p>=r[1] & p<=r[2]]
  if(length(p) > n) p <- range(p)
  }
p
}

