#' Utility Functions
#' 
#' Utility functions
#' 
#' \code{zmq.strerror()} gets ZeroMQ error message string.
#'
#' \code{zmq.version()} print current ZeroMQ version.
#' 
#' @param errno 
#' an integer for the error number
#' 
#' @return 
#' \code{zmq.strerror()} returns an R string containing ZeroMQ error
#' message.
#' 
#' @author Wei-Chen Chen \email{wccsnow@@gmail.com}.
#' 
#' @references ZeroMQ/4.1.0 API Reference:
#' \url{http://api.zeromq.org/4-1:_start}
#' 
#' Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' 
#' @examples
#' \dontrun{
#' library(pbdZMQ, quietly = TRUE)
#' 
#' context <- zmq.ctx.new()
#' zmq.ctx.destroy(context)
#' zmq.strerror(0)
#' 
#' zmq.ctx.destroy(context) # Error since context is free.
#' zmq.strerror(14)
#' }
#' 
#' @keywords internal
#' @seealso \code{\link{zmq.ctx.new}()}, \code{\link{zmq.ctx.destroy}()},
#' \code{\link{zmq.socket}()}, \code{\link{zmq.close}()}.
#' @rdname xx_utility
#' @name Utility Functions
NULL



#' @rdname xx_utility
#' @export
zmq.strerror <- function(errno){
  .Call("R_zmq_strerror", as.integer(errno[1]), PACKAGE = "pbdZMQ")
}


#' @rdname xx_utility
#' @export
zmq.version <- function(){
  ret <- .Call("R_zmq_version", PACKAGE = "pbdZMQ")
  package_version(ret)
}

