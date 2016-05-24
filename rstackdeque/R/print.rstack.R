#' @export
#' @title Print an rstack
#' @description Prints a summary of the contents of an rstack, including the number of elements and the top few.
#' @details Since only the top few elements are detailed, runs in \eqn{O(1)} time.
#' @param x the rstack to print.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @seealso \code{\link{as.list.rstack}} for converting an rstack into a list which can then be printed in full.
print.rstack <- function(x, ...) {
  s <- x
  cat(paste("An rstack with ", length(s), " elements. \n"))
  if(length(s) > 0) {
    if(length(s) > 6) {
      cat("Top of the stack:\n")
      str(as.list(head(s, 6)), comp.str = " ", no.list = T)
      cat("    ...")
    } else {
      str(as.list(head(s, length(s))), comp.str = "", no.list = T)
    }
  }
}