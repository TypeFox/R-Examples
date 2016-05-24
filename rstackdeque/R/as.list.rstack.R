#' @export
#' @title Convert an rstack to a list
#' 
#' @description Converts an rstack to a list, where the top of the stack becomes
#' the first element of the list, the second-from-top the second, and so on. 
#' 
#' @details Runs in \eqn{O(N)} time in the size of the stack, but the generated list is pre-allocated for efficiency.
#' @param x rstack to convert.
#' @param ... additional arguments passed to as.list after initial conversion to list.
#' @return a list containing the elements of the stack in top-to-bottom order.
#' @seealso \code{\link{as.data.frame.rstack}}
#' @examples
#' s <- rstack()
#' s <- insert_top(s, "a")
#' s <- insert_top(s, "b")
#' 
#' slist <- as.list(s)
#' print(slist)
as.list.rstack <- function(x, ...) {
  retlist <- vector("list", x$len)
  node <- x$head
  index <- 1
  while(!is.null(node)) {
    retlist[[index]] <- node$data
    node <- node$nextnode
    index <- index + 1
  }
  return(as.list(retlist, ...))
}