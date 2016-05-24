#' Pretty printing of Macaulay2 output.
#'
#' Pretty printing of Macaulay2 output.
#' 
#' @param x an object of class m2
#' @param ... additional parameters
#' @usage \method{print}{m2}(x, ...)
#' @return Invisible string of the printed object.
#' @export
#' @examples
#' \dontrun{
#' 
#' m2("13^1000")
#' 
#' }
#' 
print.m2 <- function(x, ...){
	
  ## argument checking and basic variable setting
  stopifnot(is.m2(x))  
  
  ## print
  sapply(as.list(x), print)
}



