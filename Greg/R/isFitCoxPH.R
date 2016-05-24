#' Functions for checking regression type
#' 
#' The \emph{isFitCoxPH} A simple check if object inherits either 
#' "coxph" or "crr" class indicating
#' that it is a survival function.
#' 
#' @param fit Regression object
#' @return \code{boolean} Returns \code{TRUE} if the object is of that type
#'  otherwise it returns \code{FALSE}.
#' 
#' @rdname isFitFn
#' @export 
#' @example inst/examples/isFitCoxPH_example.R
isFitCoxPH <- function(fit){
  if (any(c("coxph", "crr") %in% class(fit)))
    return (TRUE)
  
  return (FALSE)
}
