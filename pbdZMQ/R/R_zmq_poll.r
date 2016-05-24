#' Poll Functions
#' 
#' Poll functions
#' 
#' \code{zmq.poll()} initials ZMQ poll items given  ZMQ \code{socket}'s
#' and ZMQ poll \code{type}'s. Both \code{socket} and \code{type} are
#' in vectors of the same length, while \code{socket} contains socket pointers
#' and \code{type} contains types of poll.
#' See \code{\link{ZMQ.PO}()} for the possible values of
#' \code{type}. ZMQ defines several poll types and utilize
#' them to poll multiple sockets.
#' 
#' \code{zmq.poll.interrupt()} call \code{zmq.poll()} and raise an interrupt
#' signal if \code{ret[1] == -1} and \code{ret[2] == 4}.
#'
#' \code{zmq.poll.free()} frees ZMQ poll structure memory internally.
#'
#' \code{zmq.poll.length()} obtains total numbers of ZMQ poll items.
#'
#' \code{zmq.poll.get.revents()} obtains revent types from ZMQ poll item by
#' the input index..
#' 
#' @param socket 
#' a vector of ZMQ sockets
#' @param type 
#' a vector of socket types corresponding to \code{socket} argument
#' @param timeout
#' timeout for poll, see ZeroMQ manual for details
#' @param index
#' an index of ZMQ poll items to obtain revents
#' @param MC 
#' a message control, see \code{\link{ZMQ.MC}()} for details
#' 
#' @return \code{zmq.poll()} returns a ZMQ code and an errno,
#' see ZeroMQ manual for details, no error/warning/interrupt in this
#' \code{R} function, but some error/warning/interrupt may catch by
#' the \code{C} function \code{zmq_poll()}.
#' @return \code{zmq.poll.length()} returns the total number of poll items
#' @return \code{zmq.poll.get.revents()} returns the revent type
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
#' ### Using poll pattern.
#' ### See demo/mspoller.r for details.
#'
#' ### Run next in background or the other window.
#' SHELL> Rscript wuserver.r &
#' SHELL> Rscript taskvent.r &
#' SHELL> Rscript mspoller.r
#'
#' ### The mspoller.r has next.
#' library(pbdZMQ, quietly = TRUE)
#' 
#' ### Initial.
#' context <- zmq.ctx.new()
#' receiver <- zmq.socket(context, .pbd_env$ZMQ.ST$PULL)
#' zmq.connect(receiver, "tcp://localhost:5557")
#' subscriber <- zmq.socket(context, .pbd_env$ZMQ.ST$SUB)
#' zmq.connect(subscriber, "tcp://localhost:5556")
#' zmq.setsockopt(subscriber, .pbd_env$ZMQ.SO$SUBSCRIBE, "20993")
#' 
#' ### Process messages from both sockets.
#' cat("Press Ctrl+C or Esc to stop mspoller.\n")
#' i.rec <- 0
#' i.sub <- 0
#' while(TRUE){
#'   ### Set poller.
#'   zmq.poll(c(receiver, subscriber),
#'            c(.pbd_env$ZMQ.PO$POLLIN, .pbd_env$ZMQ.PO$POLLIN))
#' 
#'   ### Check receiver.
#'   if(bitwAnd(zmq.poll.get.revents(1), .pbd_env$ZMQ.PO$POLLIN)){
#'     ret <- zmq.recv(receiver)
#'     if(ret$len != -1){
#'       cat("task ventilator:", ret$buf, "at", i.rec, "\n")
#'       i.rec <- i.rec + 1
#'     }
#'   }
#' 
#'   ### Check subscriber.
#'   if(bitwAnd(zmq.poll.get.revents(2), .pbd_env$ZMQ.PO$POLLIN)){
#'     ret <- zmq.recv(subscriber)
#'     if(ret$len != -1){
#'       cat("weather update:", ret$buf, "at", i.sub, "\n")
#'       i.sub <- i.sub + 1
#'     }
#'   }
#' 
#'   if(i.rec >= 5 & i.sub >= 5){
#'     break
#'   }
#' 
#'   Sys.sleep(runif(1, 0.5, 1))
#' }
#' 
#' ### Finish.
#' zmq.poll.free()
#' zmq.close(receiver)
#' zmq.close(subscriber)
#' zmq.ctx.destroy(context)
#' }
#' 
#' @keywords programming
#' @seealso \code{\link{zmq.recv}()}, \code{\link{zmq.send}()}.
#' @rdname b3_poll
#' @name Poll Functions
NULL



#' @rdname b3_poll
#' @export
zmq.poll <- function(socket, type, timeout = -1L, MC = .pbd_env$ZMQ.MC){
  if(length(socket) != length(type)){
    stop("socket and type are of different length.")
  }

  type <- as.integer(type)
  if(!all(type %in% 1:7)){
    stop("type should be integers in 1 to 7.")
  }

  zmq.poll.free()

  ret <- .Call("R_zmq_poll", socket, type, as.integer(timeout),
               PACKAGE = "pbdZMQ")
  return(invisible(ret))
}


#' @rdname b3_poll
#' @export
zmq.poll.interrupt <- function(socket, type, timeout = -1L,
     MC = .pbd_env$ZMQ.MC){
  ret <- zmq.poll(socket, type, timeout, MC)

  if(ret[1] == -1 && ret[2] == 4){
    my.c <- structure(list(ret = ret), class = c("interrupt", "condition"))
    signalCondition(my.c)
  }

  return(invisible(ret))
}


#' @rdname b3_poll
#' @export
zmq.poll.free <- function(){
  ret <- .Call("R_zmq_poll_free", PACKAGE = "pbdZMQ")
  invisible(ret)
}

#' @rdname b3_poll
#' @export
zmq.poll.length <- function(){
  ret <- .Call("R_zmq_poll_length", PACKAGE = "pbdZMQ")
  invisible(ret)
}

#' @rdname b3_poll
#' @export
zmq.poll.get.revents <- function(index = 1L){
  if(index < 1){
    stop("index is a positive interger.")
  }
  ret <- .Call("R_zmq_poll_get_revents", as.integer(index - 1)[1],
               PACKAGE = "pbdZMQ")
  invisible(ret)
}
