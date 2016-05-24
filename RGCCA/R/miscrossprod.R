#' Given vectors \eqn{x} and \eqn{y} as arguments, the function miscrossprod() returns 
#' the cross-product \eqn{x^ty}. miscrossprod() handles missing data. 
#' @param x  A numeric vector.
#' @param y  A numeric vector.
#' @return \item{d.p}{The dot product between x and y: \eqn{x^ty}}
#' @title Cross product function for inputs with missing data.
#' @export miscrossprod

miscrossprod <- function(x,y) {
  d.p =  sum(drop(x)*drop(y), na.rm = TRUE)
  return(d.p)
}
