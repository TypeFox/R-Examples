#' @title One Argument
#' @description Does the function have just one argument?
#' 
#' @param f a function.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' one.arg(mean)
#' one.arg(one.arg)
#' 
#' @export
#' 
one.arg <- function(f) length(formals(f)) == 1