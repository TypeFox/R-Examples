#' @rdname newArgCheck
#' @export addWarning
#'

addWarning <- function(msg, argcheck)
{
  if (!"ArgCheck" %in% class(argcheck))
    stop("'argcheck' must be an object of class 'ArgCheck'")
  
  assign("n_warn", 
         get("n_warn", envir = argcheck) + 1, 
         envir = argcheck)
  assign("warn_msg", 
         c(get("warn_msg", envir = argcheck),
           msg), 
         envir = argcheck)
}
