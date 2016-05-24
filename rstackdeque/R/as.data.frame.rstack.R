#' @export
#' @title Convert an rstack to a data.frame
#' 
#' @description Converts the elements of an rstack into rows of a data.frame, if this is reasonable.
#' 
#' @details This function runs in \eqn{O(N)} time in the size of the rstack, and will only work if all 
#' elements of the stack have the same length() (e.g., same number of columns), and if any of the 
#' elements have names, then those names do not conflict (e.g., same column names where used). 
#' This is accomplished by a call to
#' \code{do.call("rbind", as.list.rstack(x))}, where \code{\link{as.list.rstack}} converts the rstack to a list
#' where the top element becomes the first element of the list.
#' @param x rstack to convert.
#' @param row.names passed on to \code{as.data.frame} before final conversion.
#' @param optional passed onto \code{as.data.frame} before final conversion.
#' @param ... passed onto \code{as.data.frame} before final conversion.
#' @return a data.frame with the first row the previous top of the stack.
#' @seealso \code{\link{as.list.rstack}} for conversion to a list and the generic \code{\link{as.data.frame}}.
#' @examples 
#' s <- rstack()
#' s <- insert_top(s, data.frame(names = c("Bob", "Joe"), ages = c(25, 18)))
#' s <- insert_top(s, data.frame(names = c("Mary", "Kate", "Ashley"), ages = c(27, 26, 21)))
#' print(s)
#' 
#' sd <- as.data.frame(s)
#' print(sd)
#' 
#' ## Elements may be similarly-named lists as well, representing individual rows:
#' s <- rstack()
#' s <- insert_top(s, list(name = "Bob", age = 25))
#' s <- insert_top(s, list(name = "Mary", age = 24))
#' print(s)
#' 
#' sd <- as.data.frame(s)
#' print(sd)
#' 
#' ## Building a stack in a loop, converting to a dataframe after the fact:
#' s <- rstack()
#' for(i in 1:1000) {
#'  if(runif(1,0,1) < 0.5) {
#'    s <- insert_top(s, data.frame(i = i, type = "sqrt", val = sqrt(i)))
#'  } else {
#'    s <- insert_top(s, data.frame(i = i, type = "log", val = log(i)))
#'  }
#'  if(i %% 100 == 0) {
#'    print(i/1000)
#'  }
#' }
#' print(head(as.data.frame(s)))
as.data.frame.rstack <- function(x, row.names = NULL, optional = FALSE, ...) {
  dlist <- lapply(as.list(x), as.data.frame, optional = optional, ...)
  uniquelens <- unique(lapply(dlist, length))
  if(length(uniquelens) > 1) {
    stop("cannot convert an rstack to a data.frame unless all elements have the same length()")
  }
  uniquenamesets <- unique(lapply(dlist, names))
  if(length(uniquenamesets) > 1) {
    stop("cannot convert an rstack to a data.frame when elements have contradictory names()")
  }
  return(as.data.frame(do.call(rbind, dlist), row.names, optional, ...))
}