#' New Poll Functions
#' 
#' New poll functions
#' 
#' \code{zmq.poll2()} initials ZMQ poll items given  ZMQ \code{socket}'s
#' and ZMQ poll \code{type}'s. Both \code{socket} and \code{type} are
#' in vectors of the same length, while \code{socket} contains socket pointers
#' and \code{type} contains types of poll.
#' See \code{\link{ZMQ.PO}()} for the possible values of
#' \code{type}. ZMQ defines several poll types and utilize
#' them to poll multiple sockets.
#'
#' \code{zmq.poll2.interrupt()} call \code{zmq.poll2()} and raise an interrupt
#' signal if \code{ret$pollret[1] == -1} and \code{ret$pollret[2] == 4}.
#'
#' \code{zmq.poll2.get.revents()} obtains revent types from ZMQ poll item by
#' the input index..
#' 
#' @param socket 
#' a vector of ZMQ sockets
#' @param type 
#' a vector of socket types corresponding to \code{socket} argument
#' @param timeout
#' timeout for poll, see ZeroMQ manual for details
#' @param MC 
#' a message control, see \code{\link{ZMQ.MC}()} for details
#' @param index
#' an index of ZMQ poll items to obtain revents
#' @param poller
#' a pointer of ZMQ poller
#' 
#' @return \code{zmq.poll2()} returns a ZMQ code, an errno, and a pollitem
#' pointer.
#' No error/warning/interrupt in this
#' \code{R} function, but some error/warning/interrupt may catch by
#' the \code{C} function \code{zmq_poll()}.
#' See ZeroMQ manual for details.
#' @return \code{zmq.poll.get.revents.new()} returns the revent type.
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
#' ### Multiple socket reader as in the ZeroMQ guide.
#' # SHELL> Rscript wuserver.r &
#' # SHELL> Rscript taskvent.r &
#' # SHELL> Rscript mspoller2.r
#' # SHELL> rm weather.ipc
#' 
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
#'   poller <- zmq.poll2(c(receiver, subscriber),
#'                       c(.pbd_env$ZMQ.PO$POLLIN, .pbd_env$ZMQ.PO$POLLIN))
#' 
#'   ### Check receiver.
#'   if(bitwAnd(zmq.poll2.get.revents(1, poller),
#'              .pbd_env$ZMQ.PO$POLLIN)){
#'     ret <- zmq.recv(receiver)
#'     if(ret$len != -1){
#'       cat("task ventilator:", ret$buf, "at", i.rec, "\n")
#'       i.rec <- i.rec + 1
#'     }
#'   }
#' 
#'   ### Check subscriber.
#'   if(bitwAnd(zmq.poll2.get.revents(2, poller),
#'              .pbd_env$ZMQ.PO$POLLIN)){
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
#' zmq.close(receiver)
#' zmq.close(subscriber)
#' zmq.ctx.destroy(context)
#' }
#' 
#' @keywords programming
#' @seealso \code{\link{zmq.recv}()}, \code{\link{zmq.send}()}.
#' @rdname b4_poll2
#' @name New Poll Functions
NULL



#' @rdname b4_poll2
#' @export
zmq.poll2 <- function(socket, type, timeout = -1L, MC = .pbd_env$ZMQ.MC){
  if(length(socket) != length(type)){
    stop("socket and type are of different length.")
  }

  type <- as.integer(type)
  if(!all(type %in% 1:7)){
    stop("type should be integers in 1 to 7.")
  }

  ret <- .Call("R_zmq_poll2", socket, type, as.integer(timeout),
               PACKAGE = "pbdZMQ")
  return(invisible(ret))
}

#' @rdname b4_poll2
#' @export
zmq.poll2.interrupt <- function(socket, type, timeout = -1L,
     MC = .pbd_env$ZMQ.MC){
  ret <- zmq.poll2(socket, type, timeout, MC)

  if(ret$pollret[1] == -1 && ret$pollret[2] == 4){
    my.c <- structure(list(ret = ret), class = c("interrupt", "condition"))
    signalCondition(my.c)
  }

  return(invisible(ret))
}


#' @rdname b4_poll2
#' @export
zmq.poll2.get.revents <- function(index = 1L, poller){
  if(index < 1){
    stop("index is a positive interger.")
  } else if(index > poller$pollret[3]){
    stop("index is too large.")
  }
  ret <- .Call("R_zmq_poll2_get_revents",
               as.integer(index - 1)[1], poller$pollitem,
               PACKAGE = "pbdZMQ")
  invisible(ret)
}

