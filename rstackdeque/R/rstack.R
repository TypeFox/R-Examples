#' @export
#' @title Create a new, empty rstack
#' 
#' @description An rstack is a "Last In, First Out" (LIFO) structure imagined as being organized from
#' top (last in) to bottom (first in), supporting efficient insertion into the 
#' top, removal from the top, and peeking/accessing the top element. All functions supported by rstacks are 
#' side-effect free.
#' 
#' @details Other handy functions supported by rstacks
#' include \code{as.list} and \code{as.data.frame} (the latter of which requires that
#' all elements can be appended to become rows of a data frame in a reasonable manner). Operations
#' are amortized \eqn{O(1)}.
#' 
#' The rstack class also supports \code{\link{rev}} - this operation is \eqn{O(N)}, and results in a copy. This 
#' means previous versions will retain their \eqn{O(1)} amortized nature (if we assume the cost of the 
#' reverse is charged to the newly created stack), at the cost of memory usage. However, this also means 
#' that if stacks are used in a non-persistent way, e.g. \code{s <- rev(s)}, then the garbage collector 
#' is free to clean up old versions of the data.
#' @seealso \code{\link{insert_top}} for insertion, \code{\link{without_top}} for removal, 
#' and \code{\link{peek_top}} for peeking.
#' 
#' @return an empty rstack.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' s <- rstack()
#' s <- insert_top(s, "a")
#' s <- insert_top(s, "b")
#' print(s)
#' 
#' sl <- without_top(s)
#' print(sl)
#' print(s)
#' 
#' b <- peek_top(s)
#' print(b)
rstack <- function() {
  s <- new.env(parent = emptyenv())
  s$head <- NULL
  s$tail <- NULL # for memoization
  s$len <- 0
  class(s) <- "rstack"
  return(s)
}