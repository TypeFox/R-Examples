#' Shooting Lasso
#'
#' Implementation of the Shooting Lasso (Fu, 1998) with variable dependent
#' penalization weights.
#'
#' The function implements the Shooting Lasso (Fu, 1998) with variable dependent
#' penalization. The arguments \code{XX} and \code{Xy} are optional and allow to use precalculated matrices which might improve performance.
#'
#' @param x matrix of regressor variables (\code{n} times \code{p} where
#' \code{n} denotes the number of observations and \code{p} the number of
#' regressors)
#' @param y dependent variable (vector or matrix)
#' @param lambda vector of length \code{p} of penalization parameters for each
#' regressor
#' @param control list with control parameters: \code{maxIter} maximal number
#' of iterations, \code{optTol} tolerance for parameter precision,
#' \code{zeroThreshold} threshold applied to the estimated coefficients
#' for numerical issues.
#' @param XX optional, precalculated matrix \eqn{t(X)*X}
#' @param Xy optional, precalculated matrix \eqn{t(X)*y}
#' @param beta.start start value for beta
#' @return \item{coefficients}{estimated coefficients by the Shooting Lasso
#' Algorithm} \item{coef.list}{matrix of coefficients from each iteration}
#' \item{num.it}{number of iterations run}
#' @references Fu, W. (1998). Penalized regressions: the bridge vs the lasso.
#' \emph{Journal of Computational and Graphical Software} 7, 397-416.
#' @keywords Lasso Shooting Lasso
#' @export


LassoShooting.fit <- function(x, y, lambda, control = list(maxIter = 1000, 
                                                           optTol = 10^(-5), zeroThreshold = 10^(-6)), XX = NULL, Xy = NULL, beta.start = NULL) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (is.null(XX)) 
    (XX <- crossprod(x))
  if (is.null(Xy)) 
    (Xy <- crossprod(x, y))
  # Start from the LS solution for beta if no beta.start is provided
  if (is.null(beta.start)) {
    # Ridge start beta <- MASS::ginv(XX + diag(as.vector(lambda), p) %*%
    # diag(1, p)) %*% Xy #
    # solve(XX+diag(as.vector(lambda))%*%diag(1,p))%*%Xy beta[is.nan(beta)]
    # <- 0 Zero-start beta <- rep(0,p) highest correlation start
    beta <- init_values(x, y, intercept = FALSE)$coef
  } else {
    beta <- beta.start
  }
  wp <- beta
  m <- 1
  XX2 <- XX * 2
  Xy2 <- Xy * 2
  
  while (m < control$maxIter) {
    beta_old <- beta
    for (j in 1:p) {
      # Compute the Shoot and Update the variable
      S0 <- sum(XX2[j, ] * beta) - XX2[j, j] * beta[j] - Xy2[j]
      if (sum(is.na(S0)) >= 1) {
        beta[j] <- 0
        next
      }
      
      if (S0 > lambda[j]) 
        beta[j] <- (lambda[j] - S0)/XX2[j, j]
      if (S0 < -1 * lambda[j]) 
        beta[j] <- (-1 * lambda[j] - S0)/XX2[j, j]
      if (abs(S0) <= lambda[j]) 
        beta[j] <- 0
    }
    # Update
    wp <- cbind(wp, beta)
    # Check termination
    if (sum(abs(beta - beta_old)) < control$optTol) {
      break
    }
    m <- m + 1
  }
  w <- beta
  w[abs(w) < control$zeroThreshold] <- 0
  return(list(coefficients = w, coef.list = wp, num.it = m))
}