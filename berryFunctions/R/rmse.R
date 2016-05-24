#' RMSE
#' 
#' Root Mean Squared Error
#' 
#' @details Formula used: \code{sqrt( sum((a-b)^2)/length(b) )}
#'
#' @return Numeric.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2014
#' @seealso \code{\link{rsquare}}
#' @references \url{http://en.wikipedia.org/wiki/Mean_squared_error}
#' @keywords univar
#' @export
#' @examples
#' 
#' x <- rnorm(20)
#' y <- 2*x + rnorm(20)
#' plot(x,y)
#' yp <- predict(lm(y~x))
#' plot(y, yp)
#' abline(a=0,b=1)
#' rmse(y,yp)
#' 
#' @param a Vector with values.
#' @param b Another vector of the same length.
#' @param quiet Should NA-removal warnings be suppressed? Helpful within functions. DEFAULT: FALSE
#' 
rmse <- function(
a,
b,
quiet=FALSE)
{
if(!(is.vector(a) & is.vector(b))) stop("input is not vectors")
if(length(a) != length(b)) stop("vectors not of equal length")
if(any(is.na(a)|is.na(b)))
   {
   Na <- which(is.na(a)|is.na(b))
   if(!quiet) warning(length(Na), " NAs were omitted from ", length(a), " data points.")
   a <- a[-Na] ; b <- b[-Na]
   } # end if NA
sqrt( sum((a-b)^2)/length(b) )
}
