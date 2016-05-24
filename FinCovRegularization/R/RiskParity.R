#' @title Risk Parity Portfolio
#'
#' @description
#' Computing a Risk Parity portfolio weights from the estimated
#' covariance matrix of return series.
#'
#' @param cov.mat an estimated p*p covariance matrix
#' @return a numerical vector containing the estimated portfolio weights
#' @examples
#' data(m.excess.c10sp9003)
#' assets <- m.excess.c10sp9003[,1:10]
#' RiskParity(cov(assets))
#' @importFrom stats optim
#' @export

RiskParity <- function(cov.mat) {
  RiskParity_Min <- function(w, cov.mat) {
    w <- c(w, 1 - sum(w))
    len <- length(w)
    w.mat <- matrix(rep(w, time = len), nrow = len, byrow = TRUE)
    diag.mat <- diag(w)
    res <- 2*len*as.vector(t(diag(w.mat%*%cov.mat%*%diag.mat))%*%diag(w.mat%*%cov.mat%*%diag.mat))-2*(sum(diag(w.mat%*%cov.mat%*%diag.mat)))^2
    return(res)
  }
  ncol <- ncol(cov.mat)
  weights <- optim(par = rep(0, ncol-1), fn = RiskParity_Min, cov.mat = cov.mat)$par
  weights <- c(weights, 1 - sum(weights))
  names(weights) <- colnames(cov.mat)
  return(round(weights,4))
}
