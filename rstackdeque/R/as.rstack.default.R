#' @export
#' @title Default method for converting to an rstack
#' @description Default method for converting to an rstack.
#' @details Elements from the input (of any type) are first converted to a list with \code{\link{as.list}}, 
#' after this an rstack of the appropriate size is created holding the elements. The element at \code{x[[1]]}
#' becomes the top of the stack.
#' @param x the input to convert to an rstack.
#' @param ... arguments to be passed to or from other methods (ignored).
#' @return a filled \code{\link{rstack}}.
#' @seealso \code{\link{rstack}} for info about rstacks, \code{\link{as.rstack}} for the generic.
as.rstack.default <- function(x, ...) {
  input <- x
  lastnode <- NULL
  listin <- rev(as.list(input))
  for(el in listin) {
    newnode <- rstacknode(el)
    newnode$nextnode <- lastnode
    lastnode <- newnode
  }
  newstack <- rstack()
  newstack$head <- lastnode
  newstack$len <- length(listin)
  return(newstack)
}