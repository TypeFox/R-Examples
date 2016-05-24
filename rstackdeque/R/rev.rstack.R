#' @export
#' @title Reverse an rstack
#' 
#' @description Returns a reversed version of an rstack, where the old last element (generally
#' inaccessible) is now the top (and thus now accessible).
#' 
#' @details This method runs in  \eqn{O(N)} in the size of the rstack, though it works behind-the-scenes 
#' for efficiency by converting the input stack
#' to a list, reversing the list, and building the result as a new rstack. The original is thus
#' left alone, preserving \eqn{O(1)} amortized time for the original (assuming the "cost" of reversing
#' is charged to the newly created stack) at the cost of additional memory usage. But, 
#' if the stack is not being used in a preserved manner, e.g. \code{s <- rev(s)}, the garbage collector
#' will be free to clean up the original data if it is no longer usable.
#' @param x rstack to reverse.
#' @return a reversed version of the rstack.
#' @seealso \code{\link{as.list.rstack}} for converting an rstack to a list.
#' @examples
#' s <- rstack()
#' s <- insert_top(s, "a")
#' s <- insert_top(s, "b")
#' s <- insert_top(s, "c")
#' 
#' r <- rev(s)
#' print(r)
#' print(s)
rev.rstack <- function(x) {
  saslist <- as.list(x)
  lastnode <- NULL
  for(el in saslist) {
    newnode <- rstacknode(el)
    newnode$nextnode <- lastnode
    lastnode <- newnode
  }
  newstack <- rstack()
  newstack$head <- lastnode
  newstack$len <- length(saslist)
  return(newstack)
}

