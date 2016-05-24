#' @export
#' @title Check if an rpqueue is empty
#' @description Check if an rpqueue is empty.
#' @details Runs in \eqn{O(1)} time.
#' @param x the rpqueue to check.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return logical vector of length 1.
#' @seealso \code{\link{empty}} for the generic function that can be used on rstacks, rdeques, and rpqueues.
empty.rpqueue <- function(x, ...) {
  if(length(x) > 0) {return(FALSE)}
  return(TRUE)
}

