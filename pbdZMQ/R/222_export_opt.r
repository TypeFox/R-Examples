#' Initial controls in pbdZMQ
#' 
#' Initial control functions
#' 
#' \code{.zmqopt_init()} initials default ZMQ controls.
#' \code{.zmqopt_get()} gets a ZMQ control.
#' \code{.zmqopt_set()} sets a ZMQ control.
#' 
#' @param envir 
#' an environment where ZMQ controls locate
#' @param val
#' a value to be set
#' @param main
#' a variable to be get from or set to
#' @param sub
#' a subvariable to be get from or set to
#' 
#' @return 
#' \code{.zmqopt_init()} initial the ZMQ control
#' at \code{envir}.
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
#' library(pbdZMQ, quietly = TRUE)
#' 
#' ls(.pbd_env)
#' rm(.pbd_env)
#' .zmqopt_init()
#' ls(.pbd_env)
#'
#' .pbd_env$ZMQ.SR$BLOCK
#' pbd_opt(bytext = "ZMQ.SR$BLOCK = 0L")
#' }
#' 
#' @keywords programming
#' @seealso \code{\link{.pbd_env}}.
#' @rdname a0_c_options
#' @name Initial Control Functions


### Get ZMQ options.
#' @export
#' @rdname a0_c_options
.zmqopt_get <- function(main, sub = NULL, envir = .GlobalEnv){
  if(!is.null(sub)){
    envir$.pbd_env[[main]][[sub]]
  } else{
    envir$.pbd_env[[main]]
  }
} # End of .zmqopt_get().

### Set ZMQ options.
#' @export
#' @rdname a0_c_options
.zmqopt_set <- function(val, main, sub = NULL, envir = .GlobalEnv){
  if(!is.null(sub)){
    envir$.pbd_env[[main]][[sub]] <- val
  } else{
    envir$.pbd_env[[main]] <- val
  }
  invisible()
} # End of .zmqopt_set().

### Initial ZMQ options.
#' @export
#' @rdname a0_c_options
.zmqopt_init <- function(envir = .GlobalEnv){
  if(!exists(".pbd_env", envir = envir)){
    envir$.pbd_env <- new.env()
  } 
  envir$.pbd_env$ZMQ.MC <- ZMQ.MC()
  envir$.pbd_env$ZMQ.SR <- ZMQ.SR()
  envir$.pbd_env$ZMQ.SO <- ZMQ.SO()
  envir$.pbd_env$ZMQ.ST <- ZMQ.ST()
  envir$.pbd_env$ZMQ.PO <- ZMQ.PO()

  invisible()
} # End of .zmqopt_init().

