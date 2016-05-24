#' @name GINV
#' @title Title.
#' @description Description.
#' @param M blabla
#' @return blabla
#' @author Willem Kruijer \email{willem.kruijer@@wur.nl}
#' @details blabla.
#' @examples blabla.
#' @export
GINV    <- function(M) {
  svdM    <- svd(M)
  return(svdM$v %*%diag(1/svdM$d)%*% t(svdM$u))
}
