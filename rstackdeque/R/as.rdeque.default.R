#' @export
#' @title Default method for converting to an rdeque
#' @description Default method for converting to an rdeque.
#' @details Elements from the input (of any type) are first converted to a list with \code{\link{as.list}}, 
#' after this an rdeque of the appropriate size is created holding the elements. The element at \code{x[[1]]}
#' becomes the front of the rdeque. Runs in time \eqn{O(n)}, in the size of the number of elements contained in the resulting rdeque.
#' @param x the input to convert to an rdeque
#' @param ... arguments to be passed to or from other methods (ignored).
#' @return a filled \code{\link{rdeque}}.
#' @seealso \code{\link{rdeque}} for info about rdeques, \code{\link{as.rdeque}} for the generic function.
as.rdeque.default <- function(x, ...) {
  input <- x
  newd <- rdeque()
  if(length(input) == 0) {
    return(newd)
  } else if(length(input) == 1) {
    newd$r <- insert_top(newd$r, as.list(input)[[1]])
    return(newd)
  } else {
    alllist <- as.list(input)
    mid <- as.integer(length(alllist)/2)
    left <- alllist[1:mid]
    right <- rev(alllist[(mid+1):length(alllist)])
    
    newd <- rdeque()
    newd$l <- as.rstack(left)
    newd$r <- as.rstack(right)
    return(newd)
  }
}