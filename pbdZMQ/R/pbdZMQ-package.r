#' Programming with Big Data -- Interface to ZeroMQ
#' 
#' ZeroMQ is a well-known library for high-performance
#' asynchronous messaging in scalable, distributed applications.  This
#' package provides high level R wrapper functions to easily utilize
#' ZeroMQ. We mainly focus on interactive client/server programming 
#' frameworks. For convenience, a minimal ZeroMQ library (4.1.0 rc1)
#' is shipped with pbdZMQ, which can be used if no system installation
#' of ZeroMQ is available.  A few wrapper functions compatible with
#' rzmq are also provided.
#' 
#' \tabular{ll}{ 
#'   Package: \tab pbdZMQ \cr
#'   Type: \tab Package \cr
#'   License: \tab GPL-3 2.0 \cr
#'   LazyLoad: \tab yes \cr
#' }
#' 
#' The install command using default \pkg{pbdZMQ}'s internal ZeroMQ library is
#' \cr \cr 
#' \code{> R CMD INSTALL pbdZMQ_0.1-0.tar.gz} \cr 
#' \code{--configure-args="--enable-internal-zmq"} 
#' \cr \cr 
#' Other available variables include 
#' \tabular{ll}{
#'   Variable \tab Default \cr 
#'   \code{ZMQ_INCLUDE} \tab \code{-I./zmqsrc/include} \cr 
#'   \code{ZMQ_LDFLAGS} \tab \code{-L./ -lzmq} \cr
#'   \code{ZMQ_POLLER} \tab \code{select} \cr
#' } 
#' See the package source file \code{pbdZMQ/configure.ac} for details.
#' 
#' For installation using an external ZeroMQ library, see the package source
#' file \code{pbdZMQ/INSTALL} for details.
#' 
#' @name pbdZMQ-package
#' @aliases pbdZMQ-package pbdZMQ
#' @docType package
#' @author Wei-Chen Chen \email{wccsnow@@gmail.com}.
#' @seealso \code{\link{zmq.ctx.new}()}, \code{\link{zmq.socket}()}.
#' 
#' @references ZeroMQ/4.1.0 API Reference:
#' \url{http://api.zeromq.org/4-1:_start}
#' 
#' Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' 
#' @keywords package
#' 
#' @rdname a0_a_pbdZMQ-package
NULL

