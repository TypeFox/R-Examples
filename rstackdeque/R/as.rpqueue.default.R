#' @export
#' @title Default method for converting to an rpqueue
#' @description Default method for converting to an rpqueue.
#' @details Elements from the input (of any type) are first converted to a list with \code{\link{as.list}}, 
#' after this an rpqueue of the appropriate size is created holding the elements. The element at \code{x[[1]]}
#' becomes the front of the rpqueue. Runs in time \eqn{O(n)}.
#' @param x the input to convert to an rpqueue.
#' @param ... arguments to be passed to or from other methods (ignored).
#' @return a filled \code{\link{rpqueue}}.
#' @seealso \code{\link{rpqueue}} for info about rpqueues, \code{\link{as.rpqueue}} for the generic function.
as.rpqueue.default <- function(x, ...) {
  retq <- rpqueue()
  for(el in as.list(x)) {
    retq <- insert_back(retq, el)
  }
  return(retq)  
}

