#' Send Receive Functions
#' 
#' Send and receive functions
#' 
#' \code{zmq.send()} is a high level R function calling ZMQ C API
#' \code{zmq_send()} sending \code{buf} out.
#'
#' \code{zmq.recv()} is a high level R function calling ZMQ C API
#' \code{zmq_recv()} receiving buffers of length \code{len} according to the
#' \code{buf.type}.
#'
#' \code{flags} see \code{\link{ZMQ.SR}()} for detail options of send and
#' receive functions.
#' 
#' \code{buf.type} currently supports \code{char} and \code{raw} which are both
#' in R object format.
#' 
#' @param socket 
#' a ZMQ socket
#' @param buf 
#' a buffer to be sent
#' @param len 
#' a length of buffer to be received, default 1024 bytes
#' @param flags 
#' a flag for the method using by \code{zmq_send} and
#' \code{zmq_recv}
#' @param buf.type 
#' buffer type to be received
#' 
#' @return \code{zmq.send()} returns number of bytes (invisible) in the sent
#' message if successful, otherwise returns -1 (invisible) and sets
#' \code{errno} to the error value, see ZeroMQ manual for details.
#' 
#' \code{zmq.recv()} returns a list (\code{ret}) containing the received buffer
#' \code{ret$buf} and the length of received buffer (\code{ret$len} which is
#' less or equal to the input \code{len}) if successful, otherwise returns -1
#' and sets \code{errno} to the error value, see ZeroMQ manual for details.
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
#' for(i.res in 1:5){
#'   buf <- zmq.recv(responder, 10L)
#'   cat(buf$buf, "\n")
#'   Sys.sleep(0.5)
#'   zmq.send(responder, "World")
#' }
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
#' for(i.req in 1:5){
#'   cat("Sending Hello ", i.req, "\n")
#'   zmq.send(requester, "Hello")
#'   buf <- zmq.recv(requester, 10L)
#'   cat("Received World ", i.req, "\n")
#' }
#' zmq.close(requester)
#' zmq.ctx.destroy(context)
#' }
#' 
#' @keywords programming
#' @seealso \code{\link{zmq.msg.send}()}, \code{\link{zmq.msg.recv}()}.
#' @rdname b0_sendrecv
#' @name Send Receive Functions
NULL



#' @rdname b0_sendrecv
#' @export
zmq.send <- function(socket, buf, flags = .pbd_env$ZMQ.SR$BLOCK){
  if(is.character(buf)){
    ret <- zmq.send.char(socket, buf[1], nchar(buf[1]), flags = flags)
  } else if(is.raw(buf)){
    ret <- zmq.send.raw(socket, buf[1], length(buf[1]), flags = flags)
  } else{
    stop("buf type should be char or raw.")
  }
  invisible(ret)
}



zmq.send.char <- function(socket, buf, len, flags = .pbd_env$ZMQ.SR$BLOCK){
  ret <- .Call("R_zmq_send_char", socket, buf, as.integer(len),
               as.integer(flags), PACKAGE = "pbdZMQ")
  invisible(ret)
}



zmq.send.raw <- function(socket, buf, len, flags = .pbd_env$ZMQ.SR$BLOCK){
  ret <- .Call("R_zmq_send_raw", socket, buf, as.integer(len),
               as.integer(flags), PACKAGE = "pbdZMQ")
  invisible(ret)
}



#' @rdname b0_sendrecv
#' @export
zmq.recv <- function(socket, len = 1024L, flags = .pbd_env$ZMQ.SR$BLOCK,
    buf.type = c("char", "raw")){
  if(buf.type[1] == "char"){
    ret <- zmq.recv.char(socket, len, flags = flags)
  } else if(buf.type[1] == "raw"){
    ret <- zmq.recv.raw(socket, len, flags = flags)
  } else{
    stop("buf type should be char or raw.")
  }
  invisible(ret)
}



zmq.recv.char <- function(socket, len, flags = .pbd_env$ZMQ.SR$BLOCK){
  ret <- .Call("R_zmq_recv_char", socket, as.integer(len), as.integer(flags),
               PACKAGE = "pbdZMQ")
  invisible(ret)
}



zmq.recv.raw <- function(socket, len, flags = .pbd_env$ZMQ.SR$BLOCK){
  ret <- .Call("R_zmq_recv_raw", socket, as.integer(len), as.integer(flags),
               PACKAGE = "pbdZMQ")
  invisible(ret)
}

