##' Index a colored vector
##'
##' @param x A colored factor
##' @param ... Unused
##' @return A colored factor
##' @method [ colored
##' @S3method [ colored
##' @author David C. Norris
##' @rdname Extract.colored
##' @seealso \code{\link{colored}}
`[.colored` <- function(x, ...){
  y <- NextMethod("[")
  attr(y,'colors') <- attr(x,'colors')
  y
}
