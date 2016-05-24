#' @rdname newArgCheck
#' @export addMessage
#'

addMessage <- function(msg, argcheck)
{
  if (!"ArgCheck" %in% class(argcheck))
    stop("'argcheck' must be an object of class 'ArgCheck'")
  
  assign("n_message", 
         get("n_message", envir = argcheck) + 1, 
         envir = argcheck)
  assign("message_msg", 
         c(get("message_msg", envir = argcheck),
           msg), 
         envir = argcheck)
}
