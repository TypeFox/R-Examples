##' Returns the vector of colors associated with the factor levels.
##' 
##' @param x A colored factor
##' @param \dots Additional arguments (unused)
##' @return A (character) vector of colors, having the same length as x
##' @method color colored
##' @S3method color colored
##' @author David C. Norris
##' @seealso \code{\link{colored}}
##' @keywords category color
##' 
color.colored <- function(x, ...){
  unlist(attr(x,'colors')[x])
}
