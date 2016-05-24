#' @export
#' @title Print an rdeque
#' @description Prints a summary of the contents of an rdeque, including the number of elements and the first and last few.
#' @details Depending on the internal state of the rdeque, this method is not gauranteed to run in \eqn{O(1)} time.
#' @param x the rdeque to print.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @seealso \code{\link{as.list.rdeque}} for converting an rdeque into a list which can then be printed in full.
print.rdeque <- function(x, ...) {
  d <- x
  cat(paste("A deque with ", length(d), " elements.\n"))
  if(length(d) > 0) {
    cat(paste("Front to back: \n"))
  }
  if(length(d$l) > 0) {
    str(as.list(head(d$l, 6)), comp.str = "$", no.list = T)
  }
  if(length(d) > 12) {
    cat("    ...\n")
  }
  if(length(d$r) > 0) {
    str(rev(as.list(head(d$r, 6))), comp.str = "$", no.list = T)
  }
}

