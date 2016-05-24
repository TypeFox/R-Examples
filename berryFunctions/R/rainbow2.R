#' Rainbow from blue to red
#' 
#' Reversed \code{\link{rainbow}} with different defaults, resulting in a color vector from blue (good) to red (bad)
#' 
#' @return A character vector of color names.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014
#' @seealso \code{\link{seqPal}} for a better palette, \code{\link{rainbow}}
#' @keywords color dplot
#' @export
#' @examples
#' 
#' plot(1:10, pch=16, cex=2, col=rainbow2(10))
#' 
#' @param n number of colors. DEFAULT: 10
#' @param s,v saturation and value as in \code{\link{rainbow}}. DEFAULT: 1
#' @param start start color. DEFAULT: 0
#' @param end end color. DEFAULT: 0.7
#' @param alpha transparency. DEFAULT: 1)
#' 
rainbow2 <- function(
n=10,
s=1, v=1,
start=0,
end=0.7,
alpha=1)
{
rev(rainbow(n=n, s=s, v=v, start=start, end=end, alpha=alpha))
}
