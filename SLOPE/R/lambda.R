#' Lambda sequences for SLOPE
#'
#' Computes \eqn{\lambda} sequences for SLOPE according to several pre-defined methods.
#'
#' @param n number of observations
#' @param p number of variables
#' @param fdr target False Discovery Rate (FDR)
#' @param method method to use for computing \eqn{\lambda} (see Details)
#'
#' @details The following methods for computing \eqn{\lambda} are supported:
#' \itemize{
#'  \item \code{bhq}: Computes sequence inspired by Benjamini-Hochberg (BHq)
#'    procedure
#'  \item \code{gaussian}: Computes modified BHq sequence inspired by
#'    Gaussian designs
#' }
#'
#' @rdname lambda
#' @export
create_lambda <- function(n, p, fdr=0.20, method=c('bhq','gaussian')) {
  impl = switch(match.arg(method),
                bhq = create_lambda_bhq,
                gaussian = create_lambda_gaussian_truncated)
  impl(n, p, fdr)
}

create_lambda_bhq <- function(n, p, fdr) {
  q = (1:p) * fdr / (2*p)
  qnorm(1 - q)
}

create_lambda_gaussian <- function(n, p, fdr) {
  w <- function(k) 1 / max(1, n - k - 1)
  lambda.bhq = create_lambda_bhq(n, p, fdr)
  lambda = rep(0,p)
  lambda[1] <- lambda.bhq[1]
  if (p >= 2) {
    sum_sq <- 0
    for (i in 2:p) {
      sum_sq <- sum_sq + lambda[i-1]^2
      lambda[i] <- lambda.bhq[i] * sqrt(1 + w(i-1) * sum_sq)
    }
  }
  return(lambda)
}

create_lambda_gaussian_truncated <- function(n, p, fdr) {
  lambda = create_lambda_gaussian(n, p, fdr)
  k = which.min(lambda)
  lambda[k:p] <- lambda[k]
  return(lambda)
}