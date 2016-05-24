#' File Transfer Functions
#' 
#' High level functions calling \code{zmq_send()} and \code{zmq_recv()}
#' to transfer a file in 200 KiB chunks.
#' 
#' @details
#' \code{zmq.sendfile()} binds a \code{ZMQ_PUSH} socket, and 
#' \code{zmq.recvfile()} connects to this with a \code{ZMQ_PULL} socket.
#' 
#' @param port 
#' A valid tcp port.
#' @param endpoint
#' A ZMQ socket endpoint.
#' @param filename
#' The name (as a string) of the in/out files.
#' @param verbose
#' logical; determines if a progress bar should be shown.
#' @param flags
#' a flag for the method used by \code{zmq_sendfile} and
#' \code{zmq_recvfile}
#' 
#' 
#' @return \code{zmq.sendfile()} and \code{zmq.recvfile()} return
#' number of bytes (invisible) in the sent message if successful,
#' otherwise returns -1 (invisible) and sets \code{errno} to the error
#' value, see ZeroMQ manual for details.
#' 
#' @author Drew Schmidt and Christian Heckendorf
#' 
#' @references ZeroMQ/4.1.0 API Reference:
#' \url{http://api.zeromq.org/4-1:_start}
#' 
#' Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' 
#' @examples
#' \dontrun{
#' ### Run the sender and receiver code in separate R sessions.
#' 
#' # Receiver
#' library(pbdZMQ, quietly = TRUE)
#' zmq.recvfile(55555, "localhost", "/tmp/outfile", verbose=TRUE)
#' 
#' # Sender
#' library(pbdZMQ, quietly = TRUE)
#' zmq.sendfile(55555, "/tmp/infile", verbose=TRUE)
#' }
#' 
#' @keywords programming
#' @seealso \code{\link{zmq.msg.send}()}, \code{\link{zmq.msg.recv}()}.
#' @rdname b1_sendrecvfile
#' @name File Transfer Functions
NULL



#' @rdname b1_sendrecvfile
#' @export
zmq.sendfile <- function(port, filename, verbose=FALSE,
                         flags = .pbd_env$ZMQ.SR$BLOCK)
{
  ctx <- zmq.ctx.new()
  socket <- zmq.socket(ctx, .pbd_env$ZMQ.ST$PUSH)
  endpoint <- address("*", port)
  zmq.bind(socket, endpoint)
  
  filesize <- as.double(file.info(filename)$size)
  send.socket(socket, filesize)
  
  ret <- .Call("R_zmq_send_file", socket, filename, as.integer(verbose),
               filesize, as.integer(flags),
               PACKAGE = "pbdZMQ")
  
  zmq.close(socket)
  zmq.ctx.destroy(ctx)
  
  invisible(ret)
}



#' @rdname b1_sendrecvfile
#' @export
zmq.recvfile <- function(port, endpoint, filename, verbose=FALSE,
                         flags = .pbd_env$ZMQ.SR$BLOCK)
{
  ctx <- zmq.ctx.new()
  socket <- zmq.socket(ctx, .pbd_env$ZMQ.ST$PULL)
  endpoint <- address(endpoint, port)
  zmq.connect(socket, endpoint)
  
  filesize <- receive.socket(socket)
  
  ret <- .Call("R_zmq_recv_file", socket, filename, as.integer(verbose),
               filesize, as.integer(flags),
               PACKAGE = "pbdZMQ")
  
  invisible(ret)
}
