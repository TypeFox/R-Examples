#' Knockoff filter lasso statistics
#' 
#' Computes the difference statistic
#'   \deqn{W_j = Z_j - \tilde{Z}_j}
#' or the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the lasso model.
#' 
#' @param X original design matrix
#' @param X_ko knockoff matrix
#' @param y response vector
#' @param method either 'glmnet' or 'lars' (see Details)
#' @param ... additional arguments specific to 'glmnet' or 'lars' (see Details)
#' @return The statistic \eqn{W}
#' 
#' @details This function can use either the \code{glmnet} or the \code{lars}
#' package to compute the lasso path. The \code{lars} package computes the lasso
#' path exactly, while \code{glmnet} approximates it using a fine grid of 
#' \eqn{\lambda}'s. For large matrics, \code{glmnet} can be much faster than
#' \code{lars}. By default, \code{glmnet} is used.
#' 
#' If \code{method} is set to \code{'glmnet'}, the \code{nlambda} parameter can 
#' be used to control the granularity of the grid of \eqn{\lambda}'s. The 
#' default value of \code{nlambda} is \code{5*p}, where \code{p} is the number
#' of columns of \code{X}.
#' 
#' @rdname knockoff.stat.lasso
#' @export
knockoff.stat.lasso_difference <- function(X, X_ko, y,
                                           method=c('glmnet','lars'), ...) {
  Z = lasso_max_lambda(cbind(X, X_ko), y, method, ...)
  p = ncol(X)
  orig = 1:p
  Z[orig] - Z[orig+p]
}

#' @rdname knockoff.stat.lasso
#' @export
knockoff.stat.lasso_signed_max <- function(X, X_ko, y,
                                           method=c('glmnet','lars'), ...) {
  Z = lasso_max_lambda(cbind(X, X_ko), y, method, ...)
  p = ncol(X)
  orig = 1:p
  pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])
}

#' @rdname lasso_max_lambda
#' @keywords internal
lasso_max_lambda_lars <- function(X, y) {
  if (!requireNamespace('lars', quietly=T))
    stop('lars is not installed', call.=F)
  
  fit <- lars::lars(X, y, normalize=F, intercept=F)
  lambda <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    entry <- fit$entry[j]
    if (entry > 0) lambda[j] <- fit$lambda[entry]
  }
  return(lambda)
}

#' @rdname lasso_max_lambda
#' @keywords internal
lasso_max_lambda_glmnet <- function(X, y, nlambda=NULL) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  
  n = nrow(X); p = ncol(X)
  if (is.null(nlambda)) nlambda = 5*p
  
  lambda_max = max(abs(t(X) %*% y)) / n
  lambda_min = lambda_max / 2e3
  k = (0:(nlambda-1)) / nlambda
  lambda = lambda_max * (lambda_min/lambda_max)^k
  
  fit <- glmnet::glmnet(X, y, intercept=F, standardize=F,
                        standardize.response=F, lambda=lambda)
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  indices <- apply(fit$beta, 1, first_nonzero)
  names(indices) <- NULL
  ifelse(is.na(indices), 0, fit$lambda[indices] * n)
}

#' Maximum lambda in lasso model
#' 
#' Computes the earliest (largest) lambda's for which predictors enter the
#' lasso model.
#' 
#' @param X matrix of predictors
#' @param y response vector
#' @return vector of maximum lambda's
#' 
#' @keywords internal
lasso_max_lambda <- function(X, y, method=c('glmnet','lars'), ...) {
  fn = switch(match.arg(method), 
              glmnet = lasso_max_lambda_glmnet,
              lars = lasso_max_lambda_lars)
  fn(X, y, ...)
}