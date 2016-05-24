#' @title Covariance Matrix Estimation by Macroeconomic Factor Model
#'
#' @description
#' Estimate covariance matrix by fitting a macroeconomic factor model
#' using time series regression
#'
#' @param assets a N*p matrix of asset returns, N indicates sample size
#'   and p indicates the dimension of asset returns
#' @param factor a numerical vector of length N, or a N*q matrix of
#'   macroeconomic factor(s), q indicates the dimension of factors
#' @return an estimated p*p covariance matrix
#' @examples
#' data(m.excess.c10sp9003)
#' assets <- m.excess.c10sp9003[,1:10]
#' factor <- m.excess.c10sp9003[,11]
#' MacroFactor.Cov(assets, factor)
#' @importFrom stats cov
#' @importFrom stats var
#' @export

MacroFactor.Cov <- function(assets, factor) {
  n <- dim(as.matrix(factor))[1]
  m <- dim(as.matrix(factor))[2]
  # Multivariate OLS estimator factor realizations
  X1 <- cbind(1, factor)
  beta <- solve(crossprod(X1)) %*% t(X1) %*% assets
  # Compute residual covariance matrix
  E <- assets - X1 %*% beta
  VarE <- diag(crossprod(E) / (n-m-1))
  # Compute Macroeconomics Factor covariance matrix
  if (m == 1) {
    COV <- var(factor) * beta[-1,] %*% t(beta[-1,]) + diag(VarE)
  }else{
    COV <- beta[-1,] %*% cov(factor) %*% t(beta[-1,]) + diag(VarE)
  }
  dimnames(COV) <- list(colnames(assets), colnames(assets))
  return(COV)
}
