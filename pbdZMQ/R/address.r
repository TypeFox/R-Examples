#' Form an Address/Endpoint
#' 
#' A notationally convenient function for forming addresses/endpoints.
#' It's a simple wrapper around the \code{paste0()} function.
#' 
#' @param host
#' The host ip address or url.
#' @param port
#' A port; necessary for all transports except ipc.
#' @param transport
#' The transport protocol.  Choices are "inproc", "ipc", "tcp", and
#' "pgm"/"epgm" for local in-process (inter-thread), local 
#' inter-process, tcp, and pgm, respectively.
#' 
#' @return
#' An address, for use with pbdZMQ functions.
#' 
#' @author Drew Schmidt
#' 
#' @examples
#' address("localhost", 55555)
#' 
#' @seealso
#' \code{\link{zmq.bind}}
#' 
#' @export
address <- function(host, port, transport="tcp")
{
  transports <- c("tcp", "inproc", "ipc", "pgm", "epgm")
  transport <- tolower(transport)
  match.arg(transport, transports)
  
  if (transport == "ipc")
  {
    if (!missing(port))
      warning("Ignoring specified port for ipc transport.")
    
    port <- ""
  }
  
  paste0(transport, "://", host, ":", port)
}


