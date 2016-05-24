#' Non-linear shrinkage estimator of population eigenvalues.
#'
#' @param X A data matrix.
#' @param k (Optional) Non-negative integer less than \code{ncol(X)}. If \code{k
#'   == 0} (default), \code{X} is assumed to contain 1 class, which will be
#'   centered. If \code{k >= 1}, \code{X} is assumed to contain \code{k}
#'   classes, each of which has already been centered.
#' @param method (Optional) The optimization routine used. Choices are
#'   \code{nlminb} (default) and \code{nloptr}.
#' @param control (Optional) A list of control parameters. Must correspond to
#'   the selected optimization method. See \code{\link[stats]{nlminb}},
#'   \code{\link[nloptr]{nloptr}} for details.
#' @return A numeric vector of length \code{ncol(X)}, containing the population
#'   eigenvalue estimates, sorted in ascending order.
#' @description The population eigenvalue estimates are computed by numerically
#'   inverting the QuEST function (see references). The starting point is the
#'   linear shrinkage estimate of the population eigenvalues, computed using
#'   \code{\link{linshrink}}.
#' @section NOTE: \code{nlminb} is usually robust and accurate, but does not
#'   allow equality constraints, so, in general, the sum of the estimated
#'   population eigenvalues is not equal to the sum of the sample eigenvalues.
#'   \code{nloptr} enforces an equality constraint to preserve the trace, but is
#'   substantially slower than \code{nlminb}. The default optimizer used for
#'   \code{nloptr} is the Augmented Lagrangian method with local optimization
#'   using LBFGS. These can be modified using the control parameter.
#' @references \itemize{ \item Ledoit, O. and Wolf, M. (2015). Spectrum
#'   estimation: a unified framework for covariance matrix estimation and PCA in
#'   large dimensions. Journal of Multivariate Analysis, 139(2) \item Ledoit, O.
#'   and Wolf, M. (2016). Numerical Implementation of the QuEST function.
#'   arXiv:1601.05870 [stat.CO] }
#' @examples
#' tau_estimate(X = matrix(rnorm(1e3, mean = 5), nrow = 50, ncol = 20))
#'
#' @export
tau_estimate <- function(X, k = 0, method = "nlminb", control = list()) {
  cat("Estimating population eigenvalues...")
  rho <- new.env()
  # initial set up
  rho$n <- nrow(X)
  rho$p <- ncol(X)
  rho$k <- k
  if (rho$k == 0) {
    rho$X <- X - tcrossprod(rep(1,rho$n), colMeans(X))
    rho$k <- 1
  }  else rho$X <- X

  if (rho$n > rho$k) {
    rho$effn <- rho$n - rho$k
  } else stop("k must be strictly less than nrow(X)")

  rho$S <- crossprod(rho$X)/rho$effn
  rho$lambda <- sort(eigen(rho$S, only.values = TRUE)$values)
  rho$lambda[rho$lambda < 0] = 0
  if (rho$p > rho$effn)
    rho$lambda[1:(rho$p - rho$effn)] <- 0

  rho$lambda_ls <- linshrink(X = rho$X, k = rho$k)

  # Rest of function depends on method chosen
  if (method == "nloptr") {
    lambda_eval_f <- function(tau)
    {
      if (is.unsorted(tau)) {
        tausort <- sort(tau)
        tauorder <- order(tau)
      }
      else {
        tausort <- tau
        tauorder <- 1:rho$p
      }
      Q <- QuEST(tausort,rho$effn)
      lambda_est <- Q$lambda
      lambda_est_J <- get_lambda_J(Q)
      lambda_est_order <- lambda_est[tauorder]
      lambda_est_J_order <- lambda_est_J[tauorder,]
      return (list("objective" = mean((lambda_est_order - rho$lambda)^2) ,
                   "gradient" = 2/rho$p * crossprod((lambda_est_J_order), (lambda_est_order - rho$lambda) ) ) )
    }
    lambda_eq <- function(tau) sum(tau - rho$lambda)

    if (length(control) == 0) {
      control <- list(algorithm = "NLOPT_LD_AUGLAG",
                            xtol_rel=1e-8,
                            maxeval=2000,
                            local_opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel=1e-5))
    }
    nlopt_out <- nloptr(x0 = rho$lambda_ls,
                        eval_f = lambda_eval_f,
                        lb = rep(min(rho$lambda),rho$p),
                        ub = rep(max(rho$lambda), rho$p),
                        eval_g_eq = lambda_eq,
                        eval_jac_g_eq = function(x) t(rep(1,rho$p)),
                        opts=control)
    return (nlopt_out$solution)
  }
  else if (method == "nlminb") {
    lambda_obj <- function(tau)
    {
      if (is.unsorted(tau)) {
        tausort <- sort(tau)
        tauorder <- order(tau)
      }
      else {
        tausort <- tau
        tauorder <- 1:rho$p
      }
      Q <- QuEST(tausort,rho$effn)
      lambda_est <- Q$lambda
      lambda_est_order <- lambda_est[tauorder]
      return (mean((lambda_est_order - rho$lambda)^2) )
    }

    lambda_gr <- function(tau)
    {
      if (is.unsorted(tau)) {
        tausort <- sort(tau)
        tauorder <- order(tau)
      }
      else {
        tausort <- tau
        tauorder <- 1:rho$p
      }
      Q <- QuEST(tausort,rho$effn)
      lambda_est <- Q$lambda
      lambda_est_order <- lambda_est[tauorder]
      lambda_est_J <- get_lambda_J(Q)
      lambda_est_J_order <- lambda_est_J[tauorder,]
      return (as.numeric( 2/Q$p * crossprod((lambda_est_J_order), (lambda_est_order - rho$lambda) ) ) )
    }

    optim_out <- nlminb(start = rho$lambda_ls,
                        objective = lambda_obj,
                        gradient = lambda_gr,
                        lower = min(rho$lambda),
                        upper = max(rho$lambda),
                        control = control)
    return (optim_out$par)
  }
  else stop("method should be nloptr or nlminb")
}
