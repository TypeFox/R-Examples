nlshrink_est <- function(Q) {
  delta <- rep(0,Q$p)
  if ((Q$p - Q$n) > Q$pzw) {
    lb <- min(Q$tau) - 1/ Q$n * sum(Q$tau) - 1
    ub <- 1/Q$p * (Q$p - Q$pw[1])*Q$t[1]
    optim_fn <- function(u) (sum(Q$pw*Q$t / (Q$t - u)) - Q$n)^2
    u_Fbar0 <- optimize(f = optim_fn, interval = c(lb,ub))$minimum
    null_adj <- u_Fbar0 / (1 - Q$c)
    delta[1:(Q$p - Q$n)] <- null_adj
  }

  lim0 <- ifelse(Q$p == Q$n & all(Q$tau > 0), Q$n / sum(Q$pw / Q$t), 0)

  out <- lapply(1:Q$numint, function(i) {
    y <- rep(0, Q$nidx[i])
    Dr <- abs(1 - Q$c*Q$dis_m_LF_list[[i]][Q$F_idx[[i]]])^2
    if (Q$p == Q$n) {
      j = (Q$x_F[[i]] == 0 & Dr == 0)
      y[j] <- lim0
      y[!j] <- Q$x_F[[i]][!j] / Dr[!j]
    } else {
      y <- Q$x_F[[i]] / Dr
    }
    F_temp <- Q$F[[i]]
    b <- Q$bins[[i]]
    F_temp_b <- F_temp[b]
    integral_ydF <- diff(F_temp) * 0.5 * (y[-1] + y [-Q$nidx[i]])
    integral_ydF2 <- (Q$quant[[i]] - F_temp_b) * (y[b] +
                        0.5 * (Q$quant[[i]] - F_temp_b) * (y[b + 1] - y[b]) / (F_temp[b+1] - F_temp_b) )
    integral_ydF3 <- rowSums(tcrossprod(rep(1,Q$nquant[i]), integral_ydF) * Q$integral_indic[[i]]) + integral_ydF2
    return (diff(integral_ydF3) * Q$p)
  })
  for (i in 1:Q$numint)
    delta[round(Q$F[[i]][1]*Q$p + 1):round(Q$F[[i]][Q$nidx[i]] * Q$p)] <- out[[i]]
  if (Q$p == Q$n)
    delta[delta < lim0] <- lim0

  return (delta)
}

#' Non-linear shrinkage estimator of population covariance matrix.
#'
#' @param X A data matrix.
#' @param k (Optional) Non-negative integer less than \code{ncol(X)}. If \code{k
#'   == 0} (default), \code{X} is assumed to contain 1 class, which will be
#'   centered. If \code{k >= 1}, \code{X} is assumed to contain \code{k}
#'   classes, each of which has already been centered.
#' @param method (Optional) The optimization routine called in
#'   \code{\link{tau_estimate}}. Choices are \code{nlminb} (default) and
#'   \code{nloptr}.
#' @param control (Optional) A list of control parameters. Must correspond to
#'   the selected optimization method. See \code{\link[stats]{nlminb}},
#'   \code{\link[nloptr]{nloptr}} for details.
#' @return A numeric positive semi-definite matrix of dimension \code{ncol(X)}.
#' @description \code{nlshrink_cov} calls \code{\link{tau_estimate}} to estimate
#'   the population eigenvalues. Note that the eigenvalues of the estimated
#'   population covariance matrix are not the same as the non-linear shrinkage
#'   estimates of the population eigenvalues. Theoretical and implementation
#'   details in references.
#' @references \itemize{ \item Ledoit, O. and Wolf, M. (2015). Spectrum
#'   estimation: a unified framework for covariance matrix estimation and PCA in
#'   large dimensions. Journal of Multivariate Analysis, 139(2) \item Ledoit, O.
#'   and Wolf, M. (2016). Numerical Implementation of the QuEST function.
#'   arXiv:1601.05870 [stat.CO] }
#' @examples
#' # generate matrix of uniform random variates
#' X <- matrix(sapply(1:20, function(b) runif(50, max=b)), nrow = 50, ncol = 20)
#' Sigma <- diag((1:20)^2/12) # true population covariance matrix
#' nlshrink_X <- nlshrink_cov(X, k=0) # compute non-linear shrinkage estimate
#' linshrink_X <- linshrink_cov(X, k=0) # compute linear shrinkage estimate
#' S <- cov(X) # sample covariance matrix
#'
#' # compare accuracy of estimators (sum of squared elementwise Euclidean distance)
#' sum((S-Sigma)^2)
#' sum((nlshrink_X - Sigma)^2)
#' sum((linshrink_X - Sigma)^2)
#'
#' @export
nlshrink_cov <- function(X, k=0, method = "nlminb", control = list()) {
  n <- nrow(X); p <- ncol(X)
  if (k == 0) {
    X <- X - tcrossprod(rep(1,n), colMeans(X))
    k = 1
  }

  if (n > k) {
    effn <- n - k
  } else stop("k must be strictly less than nrow(X)")
  S <- crossprod(X)/effn
  S_eigen <- eigen(S)
  U <- S_eigen$vectors
  lambda <- S_eigen$values
  lambdasort <- sort(lambda)
  lambdaorder <- order(lambda)
  tau_est <- tau_estimate(X = X, k = k, method = method, control = control)
  nlshrink_tau <- nlshrink_est(QuEST(tau_est, effn))
  return (U %*% (diag(nlshrink_tau[lambdaorder]) %*% t(U)) )
}
