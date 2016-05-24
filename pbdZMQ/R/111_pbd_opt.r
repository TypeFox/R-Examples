#' Set controls in pbdZMQ
#' 
#' Set control functions
#' 
#' \code{pbd_opt()} sets pbd options for ZMQ controls.
#'
#' \code{...} allows multiple options in
#' \code{envir$.pbd_env}, but only in a simple way.
#'
#' \code{bytext} allows to assign options by text in
#' \code{envir$.pbd_env}, but can assign advanced objects. For example,
#' \code{"option$suboption <- value"} will set
#' \code{envir$.pbd_env$option$suboption <- value}.
#' 
#' @param ... 
#' in argument format \code{option = value} to set
#' \code{.pbd_env$option <- value} inside the \code{envir}
#' @param bytext 
#' in text format \code{"option = value"} to set
#' \code{.pbd_env$option <- value} inside the \code{envir}.
#' @param envir 
#' by default the global environment is used.
#' 
#' @return 
#'  No value is returned.
#' 
#' @author Wei-Chen Chen \email{wccsnow@@gmail.com} and Drew Schmidt.
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
#' pbd_opt(bytext = "ZMQ.SR$BLOCK <- 0L")
#' }
#' 
#' @keywords programming
#' @seealso \code{\link{.pbd_env}}.
#' @rdname a0_b_pbd_opt
#' @name Set Control Functions

### Set pbd options.
#' @export
#' @rdname a0_b_pbd_opt
pbd_opt <- function(..., bytext = "", envir = .GlobalEnv){
  if(!exists(".pbd_env", envir = envir)){
    envir$.pbd_env <- new.env()
  } 

  arg <- list(...)
  if(length(arg) > 0){
    names.arg <- names(arg)
    if(is.null(names.arg) || any(names.arg == "")){
      stop("Options are all named.")
    }

    for(i.arg in 1:length(arg)){
      envir$.pbd_env[[names.arg[i.arg]]] <- arg[[i.arg]]
    }
  }

  if(bytext != ""){
    eval(parse(text = bytext), envir = envir$.pbd_env)
  }

  invisible()
} # End of pbd_opt().

