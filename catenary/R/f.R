#' Returns catenary value for given x
#' 
#' Takes catenary parameters and x and
#' return y. 
#'
#' @param x x value
#' @param c1 see \link{catenary}
#' @param c2 see \link{catenary}
#' @param lambda see \link{catenary}
#' @return y value
#' @author Jono Tuke, Matthew Roughan
#' @export
#' @note February 02 2013
#' @keywords internal
#' @examples
#' f(0,1,2,3)
f <- function(x,c1,c2,lambda){
  return(c1 * cosh((x-c2)/c1) + lambda)
}