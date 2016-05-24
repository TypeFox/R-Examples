#' Message Functions
#' 
#' Message functions
#' 
#' \code{zmq.msg.send()} sends an R message.
#' 
#' \code{zmq.msg.recv()} receives an R message.
#' 
#' @param rmsg 
#' an R message
#' @param socket 
#' a ZMQ socket
#' @param flags 
#' a flag for method of send and receive
#' @param serialize 
#' if serialize the \code{rmsg}
#' @param unserialize 
#' if unserialize the received R message
#' 
#' @return \code{zmq.msg.send()} returns 0 if successful, otherwise returns -1
#' and sets \code{errno} to \code{EFAULT}.
#' 
#' \code{zmq.msg.recv()} returns the message if successful, otherwise returns
#' -1 and sets \code{errno} to \code{EFAULT}.
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
#' ### Using request-reply pattern.
#' 
#' ### At the server, run next in background or the other window.
#' library(pbdZMQ, quietly = TRUE)
#' 
#' context <- zmq.ctx.new()
#' responder <- zmq.socket(context, .pbd_env$ZMQ.ST$REP)
#' zmq.bind(responder, "tcp://*:5555")
#' buf <- zmq.msg.recv(responder)
#' set.seed(1234)
#' ret <- rnorm(5)
#' print(ret)
#' zmq.msg.send(ret, responder)
#' zmq.close(responder)
#' zmq.ctx.destroy(context)
#' 
#' 
#' ### At a client, run next in foreground.
#' library(pbdZMQ, quietly = TRUE)
#' 
#' context <- zmq.ctx.new()
#' requester <- zmq.socket(context, .pbd_env$ZMQ.ST$REQ)
#' zmq.connect(requester, "tcp://localhost:5555")
#' zmq.msg.send(NULL, requester)
#' ret <- zmq.msg.recv(requester)
#' print(ret)
#' zmq.close(requester)
#' zmq.ctx.destroy(context)
#' }
#' 
#' @keywords programming
#' @seealso \code{\link{zmq.send}()}, \code{\link{zmq.recv}()}.
#' @rdname a2_message
#' @name Message Function
NULL



zmq.msg.init <- function(){
  .Call("R_zmq_msg_init", PACKAGE = "pbdZMQ")
}



zmq.msg.close <- function(msg.t){
  .Call("R_zmq_msg_close", msg.t, PACKAGE = "pbdZMQ")
}



#' @rdname a2_message
#' @export
zmq.msg.send <- function(rmsg, socket, flags = .pbd_env$ZMQ.SR$BLOCK,
                         serialize = TRUE){
  if(serialize){
    rmsg <- serialize(rmsg, NULL)
  }
  ret <- .Call("R_zmq_msg_send", rmsg, socket, as.integer(flags), PACKAGE = "pbdZMQ")
  invisible(ret)
}



#' @rdname a2_message
#' @export
zmq.msg.recv <- function(socket, flags = .pbd_env$ZMQ.SR$BLOCK,
                         unserialize = TRUE){
  rmsg <- .Call("R_zmq_msg_recv", socket, as.integer(flags), PACKAGE = "pbdZMQ")
  if(unserialize && is.raw(rmsg)){
    rmsg <- unserialize(rmsg)
  }
  rmsg
}

