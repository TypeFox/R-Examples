#' Sets of controls in pbdZMQ.
#' 
#' These sets of controls are used to provide default values in this package.
#' 
#' The elements of \code{.pbd_env$ZMQ.ST} are default values for socket types
#' as defined in `zmq.h' including \tabular{lcl}{ Elements \tab Value \tab
#' Usage \cr \code{PAIR} \tab 0L \tab socket type PAIR \cr \code{PUB} \tab 1L
#' \tab socket type PUB \cr \code{SUB} \tab 2L \tab socket type SUB \cr
#' \code{REQ} \tab 3L \tab socket type REQ \cr \code{REP} \tab 4L \tab socket
#' type REP \cr \code{DEALER} \tab 5L \tab socket type DEALER \cr \code{ROUTER}
#' \tab 6L \tab socket type ROUTER \cr \code{PULL} \tab 7L \tab socket type
#' PULL \cr \code{PUSH} \tab 8L \tab socket type PUSH \cr \code{XPUB} \tab 9L
#' \tab socket type XPUB \cr \code{XSUB} \tab 10L \tab socket type XSUB \cr
#' \code{STREAM} \tab 11L \tab socket type STREAM }
#' 
#' The elements of \code{.pbd_env$ZMQ.SO} are default values for socket
#' options as defined in `zmq.h' including 60 different values, see
#' \code{.pbd_env$ZMQ.SO} and `zmq.h' for details.
#' 
#' The elements of \code{.pbd_env$ZMQ.SR} are default values for send/recv
#' options as defined in `zmq.h' including \tabular{lcl}{ Elements \tab Value
#' \tab Usage \cr \code{BLOCK} \tab 0L \tab send/recv option BLOCK \cr
#' \code{DONTWAIT} \tab 1L \tab send/recv option DONTWAIT \cr \code{NOBLOCK}
#' \tab 1L \tab send/recv option NOBLOCK \cr \code{SNDMORE} \tab 2L \tab
#' send/recv option SNDMORE (not supported) }
#' 
#' The elements of \code{.pbd_env$ZMQ.MC} are default values for warning and
#' stop controls in R. These are not the ZeroMQ's internal default values. They
#' are defined as \tabular{lcl}{ Elements \tab Value \tab Usage \cr
#' \code{warning.at.error} \tab TRUE \tab if warn at error \cr
#' \code{stop.at.error} \tab TRUE \tab if stop at error }
#' 
#' @name ZMQ Control Environment
#' @aliases .pbd_env
#' @docType data
#' @format Objects contain several parameters for communicators and methods.
#' @author Wei-Chen Chen \email{wccsnow@@gmail.com}.
#' @references ZeroMQ/4.1.0 API Reference:
#' \url{http://api.zeromq.org/4-1:_start}
#' 
#' Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' @keywords global variables
#' @seealso \code{\link{.zmqopt_init}()}.
#' @rdname a0_b_control
NULL

### These are fake. These only be here for reference and to fool
### ``R CMD check''.
### The real one ``in practicee and runtime'' is initialed by the
### .zmqopt_init() which is always called by .OnLoad() in "zzz.r" to avoid
### overloaded and can be really accessed by users instead of sealed by R
### after loaded.
###
### WCC: DO ``NOT'' remark ``NOR'' use the next.
# .pbd_env <- new.env()
# .pbd_env$ZMQ.MC <- ZMQ.MC()
# .pbd_env$ZMQ.SR <- ZMQ.SR()
# .pbd_env$ZMQ.SO <- ZMQ.SO()
# .pbd_env$ZMQ.ST <- ZMQ.ST()
# .pbd_env$ZMQ.PO <- ZMQ.PO()

