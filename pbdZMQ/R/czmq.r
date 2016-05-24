#' The C-like ZeroMQ Interface
#' 
#' @description
#' The basic interface to ZeroMQ that somewhat models the C interface.
#' 
#' @details
#' A list of all functions for this interface is as follows:
#' 
#' \tabular{lll}{
#' \code{zmq.bind()} \tab \code{zmq.close()} \tab \code{zmqconnect()} \cr
#' \code{zmq.ctx.destroy()} \tab \code{zmq.ctx.new()} \tab \code{zmq.msg.recv()} \cr
#' \code{zmq.msg.send()} \tab \code{zmq.recv()} \tab \code{zmq.send()} \cr
#' \code{zmq.socket()} \tab \tab \cr
#' }
#' 
#' @author Wei-Chen Chen \email{wccsnow@@gmail.com}.
#' 
#' @references ZeroMQ/4.1.0 API Reference:
#' \url{http://api.zeromq.org/4-1:_start}
#' 
#' Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' 
#' @keywords czmq
#' @rdname xx_czmq_wrapper
#' @name C-like Wrapper Functions for ZeroMQ
NULL

