##' Extract the color key from a colored factor
##' 
##' Return a colored factor's key, mapping its levels to their colors.
##' Returns a named vector, mapping the factor levels (names) to their
##' associated colors (values).  Useful for constructing plot keys.
##' 
##' @param x A colored factor
##' @param \dots Additional arguments (unused)
##' @return A color key, in the form of a named vector
##' @method key colored
##' @S3method key colored
##' @author David C. Norris
##' @seealso \code{\link{colored}}
##' @keywords category color
##' 
key.colored <- function(x, ...){
  key <- attr(x,'colors')
  names(key) <- levels(x)
  key[!is.na(key)]
}
