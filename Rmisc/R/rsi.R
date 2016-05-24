#' Run Start Indices
#' 
#' Find the starting indices of runs in a vector.
#' 
#' @param x a vector of data.
#' 
#' @return a vector of indices indicating starting points for runs
#' 
#' @export
#' 
#' @examples
#' rsi(c(0,0,0,1,2,2,3,3,3,3,3,4))
#' 
rsi <- function(x) c(1,which(diff(x)!=0)+1)