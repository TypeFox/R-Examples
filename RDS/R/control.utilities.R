#' Set the class of the control list
#' 
#' This function sets the class of the control list, with the default being the
#' name of the calling function.
#' 
#' 
#' @param myname Name of the class to set. Defaults to the name of the calling
#' function.
#' @param control Control list. Defaults to the \code{control} variable in the
#' calling function.
#' @return The control list with class set.
#' @seealso check.control.class, print.control.list
#' @keywords utilities
set.control.class <- function(myname={sc <- sys.calls(); as.character(sc[[length(sc)-1]][[1]])}, control=get("control",pos=parent.frame())){
  class(control) <- c(myname, "control.list", "list")
  control
}

# Disable partial matching in control lists.
`$.control.list` <- getElement
