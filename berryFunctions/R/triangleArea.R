#' Area of a triangle
#' 
#' calculate Area of a planar triangle
#' 
#' @return Numeric
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2011
#' @seealso \code{\link{distance}}
#' @keywords spatial
#' @export
#' @examples
#' 
#' a <- c(1,5.387965,9); b <- c(1,1,5)
#' plot(a[c(1:3,1)], b[c(1:3,1)], type="l", asp=1)#; grid()
#' 
#' triangleArea(a,b)
#' #triangleArea(a,b[1:2])
#' 
#' @param x Vector with 3 values (x coordinates of triangle corners)
#' @param y Ditto for y.
#' @param digits Number of digits the result is rounded to. DEFAULT: 3)
#' 
triangleArea <- function(
x,
y,
digits=3)
{
if( !is.vector(x) | !is.vector(y) ) stop("Input must be a vector!")
if(length(x) != 3 | length(y) !=3 ) stop("Vectors must have 3 elements.")
A <- 0.5*(x[1] * (y[2] - y[3]) + x[2] * (y[3] - y[1]) + x[3] * (y[1] - y[2]))
return(round(abs(A),digits))
}
