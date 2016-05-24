#' @export
#' @title Return a version of an rstack without the top element
#' 
#' @description Simply returns a version of the given stack without the top element. Results in an error if the structure is empty.
#' The original rstack is left alone.
#' 
#' @details Runs in \eqn{O(1)} time worst case. 
#' 
#' @param s rstack to remove elements from.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return version of the stack with the top \eqn{n} elements removed.
#' @seealso \code{\link{insert_top}} for inserting elements.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' s <- rstack()
#' s <- insert_top(s, "a")
#' s <- insert_top(s, "b")
#' s <- insert_top(s, "c")
#' 
#' s2 <- without_top(s)
#' print(s2)
#' 
#' print(s)
without_top.rstack <- function(s, ...) {
  n = 1
  if(length(s) < n) {
    stop("cannot run without_top() on an empty rstack, check with empty() first")
  } 
  ## memozied
  if(!is.null(s$tail)) {return(s$tail)}
  newstack <- rstack()
  node <- s$head
  for(i in seq(1,n)) {
    node <- node$nextnode
  }
  newstack$len <- s$len - n
  newstack$head <- node
  s$tail <- newstack
  return(newstack)
}

