#' @title Covariance Matrix Estimation by Statistical Factor Model
#' 
#' @description
#' Estimate covariance matrix by fitting a statistical factor model 
#' using principle components analysis
#'
#' @param assets a matrix of asset returns
#' @param k numbers of factors, if k = 0, 
#'   automatically estimating by Kaiser method
#' @return an estimated p*p covariance matrix
#' @examples
#' data(m.excess.c10sp9003)
#' assets <- m.excess.c10sp9003[,1:10]
#' StatFactor.Cov(assets, 3)
#' @export

StatFactor.Cov <- function(assets, k = 0) {
  # Sample Covariance
  N <- dim(assets)[1]
  assets <- scale(assets, center = TRUE, scale = FALSE)
  cov.SAM <- (N-1)^(-1) * t(assets) %*% assets
  # SVD
  decomp <- svd(assets)
  if (k == 0) {
    # Kaiser method
    eigenvalue <- decomp$d
    k <- sum(eigenvalue > 1)
  }
  # Factor Loadings
  beta <- (N-1)^(-1/2) * decomp$v[,1:k] %*% diag(decomp$d[1:k])
  # Specific Variances
  VarE <- diag(diag(cov.SAM - beta %*% t(beta)))
  # Computing Covariance Matrix estimated by Statistical Factor Model 
  COV <- beta %*% t(beta) + VarE
  dimnames(COV) <- list(colnames(assets), colnames(assets))
  return(COV)
}