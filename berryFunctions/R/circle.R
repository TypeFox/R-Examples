#' Draw circle with a given radius
#' 
#' Draws a filled circle with a certain radius (in existing plot's units) using \code{\link{polygon}} and \code{\link{sin}}
#' 
#' @note If circles look like ellipsis, use plot(... asp=1)
#' @return none. Used for drawing.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2012
#' @seealso \code{\link{symbols}}, \code{\link{polygon}}
#' @keywords aplot
#' @export
#' @examples
#' 
#' plot(1:20, type="n", asp=1)
#' circle(5,5, r=3)   # 1:1 aspect shows they're really circles and not ellipses.
#' circle(15,10, r=4, locnum=12, col=2, border=4, lwd=3)
#' 
#' # can not be vectorized:
#' x <- sample(1:20, 15) ;  y <- sample(1:20, 15) ; r <- runif(20)*3
#' circle(x,y,r, col=rgb(1,0.5,0,alpha=0.4), border=NA)
#' for(i in 1:15) circle(x[i],y[i],r[i], col=rgb(1,0.5,0,alpha=0.4), border=NA)
#' 
#' @param x x coordinate of points, numeric value of length 1
#' @param y y coordinate
#' @param r radius of the circle in units of current plot
#' @param locnum number of calculated points on the circle (more means smoother but slower). DEFAULT: 100
#' @param \dots further arguments passed to \code{\link{polygon}}, like col, border, lwd
#' 
circle <- function(
  x, 
  y, 
  r, 
  locnum=100, 
  ...) 
{
# input checking - only one circle can be drawn:
if(length(x) >1 | length(y) >1 | length(r) >1 | length(locnum) >1)
  {
  warning("Only the first element of the vectors is used.")
  x <- x[1]; y <- y[1]; r <- r[1]; locnum <- locnum[1]
  }
# input checking - is every value numeric?
if(!is.numeric(x)) stop("x must be numeric")
if(!is.numeric(y)) stop("y must be numeric")
if(!is.numeric(r)) stop("r must be numeric")
# prepare circle line coordinates:
cx <- x+r*cos( seq(0,2*pi,len=locnum) )
cy <- y+r*sin( seq(0,2*pi,len=locnum) )
# actually draw it:
polygon(cx, cy, ...)
}
