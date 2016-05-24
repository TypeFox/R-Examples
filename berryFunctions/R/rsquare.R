#' Rsquare
#' 
#' R squared (coefficient of determination)
#' 
#' @details Formula used: \code{\link{cor}(a,b)^2}
#' 
#' @return Numeric.
#' @note Using cor is much faster than using\cr \code{ aa <- a-mean(a); bb <- b-mean(b); sum(aa*bb)^2/sum(aa^2)/sum(bb^2)}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2014
#' @seealso \code{\link{rmse}}, \code{\link{cor}}, \code{\link{lm}}
#' @references \url{http://en.wikipedia.org/wiki/R-squared}
#' @keywords univar
#' @export
#' @examples
#' 
#' x <- rnorm(20)
#' y <- 2*x + rnorm(20)
#' plot(x,y)
#' rsquare(x,y)
#' 
#' r2 <- sapply(1:10000, function(i){
#'    x <- rnorm(20);  y <- 2*x + rnorm(20);  rsquare(x,y) })
#' hist(r2, breaks=70, col=5,
#' main= "10'000 times   x <- rnorm(20);  y <- 2*x + rnorm(20);  rsquare(x,y)")
#' 
#' @param a Vector with values.
#' @param b Another vector of the same length.
#' @param quiet Should NA-removal warnings be suppressed? Helpful within functions. DEFAULT: FALSE
#' 
rsquare <- function(
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
if(all(a==0) | all(b==0))
  {
  if(!quiet) warning("all a (or all b) values are zero, returning NA.")
  return(NA)
  } # end if zero
cor(a,b)^2
}



if(FALSE)
{
# alternative, slower (3.4 instead of 2.1 seconds in the example below)
# crucial, if calculations are done iteratively or performed multiple times
rsquare2 <- function(a,b) { 
  if(!(is.vector(a) & is.vector(b))) stop("input is not vectors")
  if(length(a) != length(b)) stop("vectors not of equal length")
  if(any(is.na(a)|is.na(b)))
     { warning("NAs were omitted")
     Na <- which(is.na(a)|is.na(b))
     a <- a[-Na] ; b <- b[-Na]
     } # end if NA
  aa <-  a-mean(a)
  bb <-  b-mean(b)
  sum(aa*bb)^2/sum(aa^2)/sum(bb^2) }

a <- sort(rnorm(1e8)); b <- 2*a+3+rnorm(length(a))
system.time(rsquare(a,b))
system.time(rsquare2(a,b))
}
