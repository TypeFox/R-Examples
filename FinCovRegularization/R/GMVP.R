#' @title Global Minimum Variance Portfolio
#'
#' @description
#' Computing a global minimum variance portfolio weights from the estimated 
#' covariance matrix of return series.
#' 
#' @param cov.mat an estimated p*p covariance matrix
#' @param short logical flag, indicating whether shortsales on the 
#'   risky assets are allowed
#' @return a numerical vector containing the estimated portfolio weights
#' @examples
#' data(m.excess.c10sp9003)
#' assets <- m.excess.c10sp9003[,1:10]
#' GMVP(cov(assets), short=TRUE)
#' GMVP(cov(assets), short=FALSE)
#' @export

GMVP <- function(cov.mat, short=TRUE) {
  n <- ncol(cov.mat)
  dvec <- matrix(0,n,1)
  if (short==TRUE) {
    Amat <- matrix(rep(1,n), n, 1)
    bvec <- 1
  } else if (short == FALSE) {
    Amat <- matrix(c(rep(1,n),diag(n)), n, n+1)
    bvec <- c(1, rep(0,n))
  }
  meq <- 1
  weights <- quadprog::solve.QP(Dmat = cov.mat, dvec, Amat, bvec, meq)$solution
  names(weights) <- colnames(cov.mat)
  return(round(weights, 4))
}