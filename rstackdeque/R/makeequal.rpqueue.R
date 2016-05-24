#' @export
#' @title Maintenance function for rpqueues
#' @description Maintenance function for rpqueues, called automatically when needed by other functions.
#' @details See \emph{Simple and Efficient Purely Functional Queues and Deques}, 
#' Okasaki 1995, J. Functional Programming, 5(4) 583 to 592 for information.
#' @param rpqueue rpqueue to makeequal.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return a "fixed" rpqueue.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @seealso \code{\link{rotate}} helper function that calls this one.
makeequal.rpqueue <- function(rpqueue, ...) {
  #print(paste(length(rpqueue$lhat, length(rpqueue$l), length(rpqueue$r))))
  if(length(rpqueue$lhat) > 0) {
    # maybe here we can avoid creating a whole new object?
    newq <- rpqueue()
    newq$l <- rpqueue$l
    newq$r <- rpqueue$r
    newq$lhat <- without_top(rpqueue$lhat)
    return(newq)
  } else {
    newq <- rpqueue()
    acc <- rstack()
    resq <- rotate(rpqueue, acc)
    newq$l <- resq$l
    newq$lhat <- resq$l
    newq$r <- rstack()
    return(newq)
  }
}