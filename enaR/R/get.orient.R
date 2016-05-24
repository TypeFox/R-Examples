#' get.orient --- returns the global orientation
#' INPUT = none
#' OUTPUT = returns the current orientation of matrices
#' 
#' M. Lau | 16 Jun 2013
#' ---------------------------------------------------

get.orient <- function(){
  current.orientation <- get('orientation',envir=environment(set.orient))
  return(current.orientation)
}
