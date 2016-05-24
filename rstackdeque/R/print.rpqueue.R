#' @export
#' @title Print an rpqueue
#' @description Prints a summary of the contents of an rpqueue, including the number of elements and the first few.
#' @details Since only the first few elements are detailed, runs in \eqn{O(1)} time.
#' @param x the rpqueue to print.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @seealso \code{\link{as.list.rpqueue}} for converting an rpqueue into a list which can then be printed in full.
print.rpqueue <- function(x, ...) {
  cat(paste("A queue with ", length(x), " elements.\n"))
  if(length(x) > 0) {
    cat(paste("Front: \n"))
  }
  if(length(x) > 0) {
    str(as.list(head(x, 6)), comp.str = "$", no.list = T)
  }
}