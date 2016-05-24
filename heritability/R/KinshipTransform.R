#' @name KinshipTransform
#' @title Title.
#' @description Description.
#' @param Matrix blabla
#' @return blabla
#' @author Willem Kruijer \email{willem.kruijer@@wur.nl}
#' @details blabla.
#' @examples blabla.
#' @export
KinshipTransform <- function(Matrix) {
  nn <- ncol(Matrix)
  return((sum(diag(Matrix))-as.numeric(matrix(1,1,nn) %*% Matrix %*% matrix(1,nn,1))/nn)/(nn-1))
}
