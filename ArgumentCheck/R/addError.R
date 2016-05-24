#' @rdname newArgCheck
#' @export addError
#'

addError <- function(msg, argcheck)
{
  if (!"ArgCheck" %in% class(argcheck))
    stop("'argcheck' must be an object of class 'ArgCheck'")
  
  assign("n_error", 
         get("n_error", envir = argcheck) + 1, 
         envir = argcheck)
  assign("error_msg", 
         c(get("error_msg", envir = argcheck),
           msg), 
         envir = argcheck)
}

