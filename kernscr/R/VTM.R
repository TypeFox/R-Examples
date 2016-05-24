#'VTM
#'
#'Repeats a vector as identical rows of a matrix
#'
#'@param vc vector to be repeated
#'@param dm number of times that \code{vc} should be repeated
#'
#'@return a matrix of \code{dm} rows
#'
#'@keywords internal
#'@export
VTM <- function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}