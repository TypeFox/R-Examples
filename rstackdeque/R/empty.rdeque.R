#' @export
#' @title Check if an rdeque is empty
#' @description Check if an rdeque is empty.
#' @details Runs in \eqn{O(1)} time.
#' @param x the rdeque to check.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return logical vector of length 1.
#' @seealso \code{\link{empty}} for the generic function that can be used on rstacks, rdeques, and rpqueues.
empty.rdeque <- function(x, ...) {
  if(length(x) < 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


