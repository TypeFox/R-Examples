#' Random Port
#' 
#' Generate a valid, random TCP port.
#' 
#' By definition, a TCP port is an unsigned short, and so it can not
#' exceed 65535.  Additionally, ports in the range 1024 to 49151 are
#' (possibly) registered by ICANN for specific uses.
#' 
#' \code{random_port()} will simply generate a valid, non-registered
#' tcp port.  \code{random_unused_port()} will generate a port
#' that is available for socket connections.
#' 
#' \code{random_open_port()} finds a random port not already bound
#' to an endpoint.
#' 
#' @param min_port,max_port
#' The minimum/maximum value to be generated.  The minimum should not
#' be below 49152 and the maximum should not exceed 65536 (see
#' details).
#' @param max_tries
#' The maximum number of times a random port will be searched for.
#' 
#' @references
#' "The Ephemeral Port Range" by Mike Gleason.  
#' \url{http://www.ncftp.com/ncftpd/doc/misc/ephemeral_ports.html}
#' 
#' @author Drew Schmidt
#' 
#' @examples
#' random_port()
#' 
#' @importFrom stats runif
#' 
#' @rdname random_port
#' @export
random_port <- function(min_port=49152, max_port=65536)
{
  if (min_port < 49152)
    warning("non-recommended 'min_port' value; see ?random_port for details")
  if (max_port > 65536)
    warning("non-recommended 'max_port' value; see ?random_port for details")
  
  
  as.integer(runif(1, min_port, max_port+1))
}



#' @rdname random_port
#' @export
random_open_port <- function(min_port=49152, max_port=65536, max_tries=100)
{
  ctxt <- init.context()
  socket <- init.socket(ctxt, "ZMQ_REP")
  
  ret <- -1
  
  for (i in 1:max_tries)
  {
    port <- random_port(min_port=min_port, max_port=max_port)
    addr <- paste0("tcp://*:", port)
    
    catch <- tryCatch(zmq.bind(socket, addr), error=identity, warning=identity)
    if (!inherits(catch, "error") && !inherits(catch, "warning"))
    {
      ret <- port
      break
    }
  }
  
  zmq.close(socket)
  zmq.ctx.destroy(ctxt)
  rm(socket);rm(ctxt);
  invisible(gc)
  
  if (ret == -1)
    stop("No valid port could be found.")
  else
    return(ret)
}


