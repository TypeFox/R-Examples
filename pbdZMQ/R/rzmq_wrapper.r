#' All Wrapper Functions for rzmq
#' 
#' Wrapper functions for backwards compatibility with rzmq.  See vignette
#' for examples.
#' 
#' @details
#' \code{send.socket()}/\code{receive.socket()} send/receive messages over
#' a socket.  These are simple wrappers around \code{zmq.msg.send()} and
#' \code{zmq.msg.receive()}, respectively.
#' 
#' \code{init.context()} creates a new ZeroMQ context.  A useful wrapper
#' around \code{zmq.ctx.new()} which handles freeing memory for you, i.e.
#' \code{zmq.ctx.destroy()} will automatically be called for you.
#' 
#' \code{init.socket()} creates a ZeroMQ socket; serves as a high-level
#' binding for \code{zmq.socket()}, including handling freeing memory
#' automatically.  See also \code{.pbd_env$ZMQ.ST}.
#' 
#' \code{bind.socket()}:  see \code{zmq.bind()}.
#' 
#' \code{connect.socket()}:  see \code{zmq.connect()}
#' 
#' @param socket
#' A ZMQ socket.
#' @param data
#' An R object.
#' @param send.more
#' Logical; will more messages be sent?
#' @param serialize,unserialize
#' Logical; determines if serialize/unserialize should be called
#' on the sent/received data.
#' @param dont.wait
#' Logical; determines if reception is blocking.
#' @param context
#' A ZMQ context.
#' @param socket.type
#' The type of ZMQ socket as a string, of the form "ZMQ_type".  Valid 'type'
#' values are PAIR, PUB, SUB, REQ, REP, DEALER, PULL, PUSH, XPUB, XSUB, and
#' STERAM.
#' @param address
#' A valid address.  See details.
#' 
#' @author Wei-Chen Chen \email{wccsnow@@gmail.com}.
#' 
#' @references ZeroMQ/4.1.0 API Reference:
#' \url{http://api.zeromq.org/4-1:_start}
#' 
#' Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' 
#' @keywords rzmq
#' @rdname xx_rzmq_wrapper
#' @name Wrapper Functions for rzmq
NULL



#' @rdname xx_rzmq_wrapper
#' @export
send.socket <- function(socket, data, send.more = FALSE, serialize = TRUE){
  if(send.more){
    flags <- .pbd_env$ZMQ.SR$SNDMORE
  } else{
    flags <- .pbd_env$ZMQ.SR$BLOCK
  }
  zmq.msg.send(data, socket, flags = flags, serialize = serialize)
}



#' @rdname xx_rzmq_wrapper
#' @export
receive.socket <- function(socket, unserialize = TRUE, dont.wait = FALSE){
  if(dont.wait){
    flags <- .pbd_env$ZMQ.SR$DONTWAIT
  } else{
    flags <- .pbd_env$ZMQ.SR$BLOCK
  }
  zmq.msg.recv(socket, flags = flags, unserialize = unserialize)
}



#' @rdname xx_rzmq_wrapper
#' @export
init.context <- function(){
  try.zmq.ctx.destroy <- function(ctx){
    invisible(suppressWarnings(zmq.ctx.destroy(ctx)))
  }
  
  ctx <- zmq.ctx.new()  
  reg.finalizer(ctx, try.zmq.ctx.destroy, onexit = TRUE)
  ctx 
}



#' @rdname xx_rzmq_wrapper
#' @export
init.socket <- function(context, socket.type){
  try.zmqt.close <- function(socket){
    invisible(suppressWarnings(zmq.close(socket)))
  }
  
  socket.type <- sub(".*_", "", socket.type)
  id <- which(names(.pbd_env$ZMQ.ST) == socket.type)
  if(length(id) != 1){
    stop("socket.type is not found.")
  } else{
    socket.type <- .pbd_env$ZMQ.ST[[id]]
  }
  
  socket <- zmq.socket(context, type = socket.type)
  reg.finalizer(socket, try.zmqt.close, onexit = TRUE)
  socket
}



#' @rdname xx_rzmq_wrapper
#' @export
bind.socket <- function(socket, address){
  zmq.bind(socket, address)
}



#' @rdname xx_rzmq_wrapper
#' @export
connect.socket <- function(socket, address){
  zmq.connect(socket, address)
}

