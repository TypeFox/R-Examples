#' @export
#' @title Maintenance function for rpqueues
#' @description Maintenance function for rpqueues, called automatically when needed by other functions.
#' @details See \emph{Simple and Efficient Purely Functional Queues and Deques}, 
#' Okasaki 1995, J. Functional Programming, 5(4) 583 to 592 for information on this function.
#' @param rpqueue rpqueue to rotate.
#' @param acclazylist lazy list accumulator.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return a fully rotated rpqueue, but with the l list delayed in evaluation.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @seealso \code{\link{makeequal}} helper function that calls this one.
rotate.rpqueue <- function(rpqueue, acclazylist, ...) {
  ## safety check :-P
  if(length(rpqueue) == 0) {
    return(rpqueue)
  }
  if(length(rpqueue$l) == 0) {
    newq <- rpqueue()
    newq$l <- rstack()
    
    newq$l$head <- rstacknode(peek_top(rpqueue$r))
    newq$l$head$nextnode <- acclazylist$head  #ie, tail is the lazylist
    newq$l$len <- length(rpqueue) + length(acclazylist)
    newq$r <- rstack()
    return(newq)
  } else {
    newq <- rpqueue()
    newq$l <- rstack()
    newq$l$head <- rstacknode(peek_top(rpqueue$l))
    newq$l$len <- length(rpqueue) + length(acclazylist)
    
    without_heads <- rpqueue()
    without_heads$l <- without_top(rpqueue$l)
    without_heads$r <- without_top(rpqueue$r)
    
    acc <- insert_top(acclazylist, peek_top(rpqueue$r))
    delayedAssign("nextnode", rotate(without_heads, acc)$l$head, assign.env = newq$l$head)
    newq$r <- rstack()
    return(newq)
  }
}