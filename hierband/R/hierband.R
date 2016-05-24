#' Solves main optimization problem for fixed lambda value
#'
#' Solves the main optimization problem appearing in Bien, Bunea, & Xiao (2015):
#' \deqn{min_P || Sighat - P ||_F^2 + lam * sum_l || (W_l * P)_{g_l} ||_2}
#' where g_l are the outermost l(l+1) elements of a square matrix.
#' and || (W_l * P)_{g_l} ||^2 = sum_{m<=l} w_{lm}^2 ||P_{s_m}||^2.  If a non-NULL \code{delta} is provided,
#' then a constraint of the form $P >= delta I_p$ is included. Problem is solved by
#' performing blockwise coordinate descent on the dual problem.  See paper
#' for more explanation.
#'
#' @param Sighat The sample covariance matrix
#' @param lam Non-negative penalty parameter.  Controls sparsity level.
#' @param w \code{(p-1)}-by-\code{(p-1)} lower-triangular matrix (above diagonal ignored).
#'          \code{w[l,]} gives the \code{l} weights for g_l.
#'          Defaults to \code{w[l,m]=sqrt(2 * l)/(l - m + 1)} for \code{m <= l}
#' @param delta Lower bound on eigenvalues.  If this is NULL (which is default), then no eigenvalue
#' constraint is included.
#' @param maxiter Number of iterations of blockwise coordinate descent to perform.
#' @param tol Only used when \code{delta} is non-NULL.  When no eigenvalue changes by more than
#' \code{tol} in BCD, convergence is assumed.
#'
#' @return Returns the convex banded estimate of covariance.
#'
#' @export
#'
#' @seealso \code{\link{hierband.path}} \code{\link{hierband.cv}}
#' @examples
#' set.seed(123)
#' p <- 100
#' n <- 50
#' K <- 10
#' true <- ma(p, K)
#' x <- matrix(rnorm(n*p), n, p) %*% true$A
#' Sighat <- cov(x)
#' fit <- hierband(Sighat, lam=0.4)
#' min(eigen(fit)$values)
#' fit2 <- hierband(Sighat, lam=0.4, delta=0.2)
#' min(eigen(fit2)$values)
#' # Use cross validation to select lambda:
#' path <- hierband.path(Sighat)
#' cv <- hierband.cv(path, x)
#' fit <- hierband(Sighat, lam=cv$lam.best)
hierband <- function(Sighat, lam, w=NULL, delta=NULL, maxiter=100, tol=1e-7) {
  p <- nrow(Sighat)
  stopifnot(ncol(Sighat)==p)
  stopifnot(Sighat==t(Sighat))
  if (is.null(w)) {
    #w <- matrix(NA, p-1, p-1)
    #for (l in seq(p-1)) for (m in seq(l)) w[l, m] <- sqrt(2 * l) / (l - m + 1)
    w <- formw(p)
  }
  else stopifnot(is.matrix(w), nrow(w)==p-1, ncol(w)==p-1)
  stopifnot(lam>=0, w[lower.tri(w)]>=0)
  if (is.null(delta)) {
    # no eigenvalue constraint
    return(subdiag.thresh(R=Sighat, lam=lam, w=w))
  } else {
    # eigenvalue constraint
    lamC <- matrix(0, p, p)
    for (i in seq(maxiter)) {
      R.A <- Sighat + lamC # partial residual excluding A terms
      R <- subdiag.thresh(R=R.A, lam=lam, w=w)
      eig <- eigen(lamC - R) # partial residual excluding C term
      evals.new <- pmax(eig$values + delta, 0)
      if (i > 1) {
        if (max(abs(evals.new-evals)) < tol) {
          cat("Converged after", i, "iterations.", fill=TRUE)
          break
        }
      }
      evals <- evals.new
      lamC <- eig$vec %*% diag(evals) %*% t(eig$vec)
    }
  }
  R
}


#' Performs a single pass of BCD on a matrix R.
#'
#' To solve the unconstrained problem, R is Sigmahat.
#' To solve constrained problem, R is the current partial residual (excluding A terms).
#' 
#' @param R p-by-p symmetric matrix
#' @param lam Non-negative penalty parameter.  Controls sparsity level.
#' @param w \code{(p-1)}-by-\code{(p-1)} lower-triangular matrix (above diagonal ignored).
#'          \code{w[l,]} gives the \code{l} weights for g_l.
#'          Defaults to \code{w[l,m]=sqrt(2 * l)/(l - m + 1)} for \code{m <= l}
subdiag.thresh <- function(R, lam, w=NULL) {
  SMALL <- 1e-14
  p <- nrow(R)
  if (is.null(w)) {
    #w <- matrix(NA, p-1, p-1)
    #for (l in seq(p-1)) for (m in seq(l)) w[l, m] <- sqrt(2 * l) / (l - m + 1)
    w <- formw(p)
  } else {
    if (all( w == diag(diag(w)) ))
      return(gpband(R, lam, diag(w)))
    if (any(w[lower.tri(w)] == 0))
      stop("w must either have all positive weights or else be diagonal.")
  }
  #T <- toeplitz(p:1)
  #r <- as.numeric(sqrt(tapply(R^2,T,sum)))[-p] # r[l] = ||R_{s_l}||_2
  r <- sqrt(2) * rev(subdiagonal.l2norms(R)[-1])
  # multiplied by sqrt(2) since s_\ell includes two identical subdiags
  pp <- r
  nu <- rep(NA, p-1)
  for (l in seq(p-1)) {
      ww <- w[l,1:l]
      if (sum((pp[1:l]/ww)^2) <= lam^2) {
        nu[l] <- 0
        pp[1:l] <- 0
      } else {
        f <- function(nu) 1 - lam / sqrt(sum(( pp[1:l]/(ww + nu/ww) )^2)) # according to Qin... Goldfarb, this is stabler
        nu.u <- sqrt(sum((ww*pp[1:l])^2)) / lam
        nu.l <- max(nu.u - ww[l]^2, 0)
        if (abs(f(nu.l)) < SMALL)
          nu[l] <- nu.l
        else if (abs(f(nu.u)) < SMALL)
          nu[l] <- nu.u
        else {
          nu[l] <- uniroot(f, c(nu.l, nu.u), tol=1e-15)$root
        }
        pp[1:l] <- nu[l] * pp[1:l] / (ww^2 + nu[l])
      }
  }
  That <- toeplitz(c(1,rev(pp / r))) # adaptive tapering matrix
  That * R # returns updated residual matrix
}

#' Groupwise soft-thresholds subdiagonals by lam * w
#'
#' @param R p-by-p symmetric matrix
#' @param lam Non-negative penalty parameter.  Controls sparsity level.
#' @param w \code{(p-1)} vector of weights. Default, w[l]=sqrt(2 * l)
gpband <- function(R, lam, w=NULL) {
  # Solves min_P || R - P ||_F^2 + lam * sum_l w_l || P_{l} ||_2
  # where P_l is the (two) lth diagonals
  p <- nrow(R)
  stopifnot(R==t(R))
  if (is.null(w)) w <- sqrt(2 * seq(p-1))
  stopifnot(lam >= 0, w >= 0)
  r <- sqrt(2) * subdiagonal.l2norms(R)[-1]
  tt <- pmax(1 - lam * rev(w) / r, 0)
  toeplitz(c(1, tt)) * R
}


#' Compute the L2 norm of each subdiagonal of a symmetric matrix R.
#'
#' @param R a symmetric matrix
#' @export
#' @useDynLib hierband
subdiagonal.l2norms <- function(R) {
  p <- nrow(R)
  stopifnot(ncol(R)==p)
  out <- .C("subdiag_l2norm",
            as.double(R),
            as.integer(p),
            r=rep(0, p), PACKAGE="hierband")
  out$r
}

#' Form the "general weights" matrix 
#'
#' @param p dimension of covariance matrix
#' @export
#' @useDynLib hierband
formw <- function(p) {
  out <- .C("formw",
            as.integer(p),
            w=rep(0,(p-1)^2),
            PACKAGE="hierband")
  matrix(out$w, nrow=p-1, byrow = TRUE)
}
