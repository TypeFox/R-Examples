#' Send Receive Multiple Raw Buffers
#' 
#' Send and receive functions for multiple raw buffers
#' 
#' \code{zmq.send.multipart()} is a high level R function to send multiple
#' raw messages \code{parts} at once.
#' 
#' \code{zmq.recv.multipart()} is a high level R function to receive multiple
#' raw messages at once.
#' 
#' @param socket 
#' a ZMQ socket
#' @param parts
#' a vector of multiple buffers to be sent
#' @param serialize,unserialize
#' if serialize/unserialize the received multiple buffers
#' 
#' @return \code{zmq.send.multipart()} returns.
#'
#' \code{zmq.recv.multipart()} returns.
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
#' 
#' ret <- zmq.recv.multipart(responder, unserialize = TRUE)
#' parts <- as.list(rep("World", 5))
#' zmq.send.multipart(responder, parts)
#' for(i in 1:5) cat(ret[[i]])
#' 
#' zmq.close(responder)
#' zmq.ctx.destroy(context)
#' 
#' ### At a client, run next in foreground.
#' library(pbdZMQ, quietly = TRUE)
#' 
#' context <- zmq.ctx.new()
#' requester <- zmq.socket(context, .pbd_env$ZMQ.ST$REQ)
#' zmq.connect(requester, "tcp://localhost:5555")
#' 
#' parts <- lapply(1:5, function(i.req){ paste("Sending Hello ", i.req, "\n") })
#' zmq.send.multipart(requester, parts)
#' ret <- zmq.recv.multipart(requester, unserialize = TRUE)
#' print(ret)
#' 
#' zmq.close(requester)
#' zmq.ctx.destroy(context)
#' }
#' 
#' @keywords programming
#' @seealso \code{\link{zmq.msg.send}()}, \code{\link{zmq.msg.recv}()}.
#' @rdname b2_sendrecv_multipart
#' @name Send Receive Multiple Raw Buffers
NULL



#' @rdname b2_sendrecv_multipart
#' @export
zmq.send.multipart <- function(socket, parts, serialize = TRUE){
  for(i.part in 1:(length(parts) - 1)){
    zmq.msg.send(parts[[i.part]], socket, flags = .pbd_env$ZMQ.SR$SNDMORE,
                 serialize = serialize)
  }
  zmq.msg.send(parts[[length(parts)]], socket, flags = .pbd_env$ZMQ.SR$BLOCK,
               serialize = serialize)

  invisible()
}



#' @rdname b2_sendrecv_multipart
#' @export
zmq.recv.multipart <- function(socket, unserialize = TRUE){
  ret <- list() 
  i.part <- 1
  ret[[i.part]] <- zmq.msg.recv(socket, flags = .pbd_env$ZMQ.SR$BLOCK,
                                unserialize = unserialize)
  opt.val <- zmq.getsockopt(socket, .pbd_env$ZMQ.SO$RCVMORE, 0L)

  while(opt.val == 1){
    i.part <- i.part + 1
    ret[[i.part]] <- zmq.msg.recv(socket, flags = .pbd_env$ZMQ.SR$BLOCK,
                                  unserialize = unserialize)
    opt.val <- zmq.getsockopt(socket, .pbd_env$ZMQ.SO$RCVMORE, 0L)
  }

  ret
}

