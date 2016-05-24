#' Solves main optimization problem over a grid of lambda values
#' 
#' See \code{\link{hierband}} for the problem this is solving.  If \code{lamlist} not provided, then grid will be constructed
#' starting at lambda_max, the smallest value of lam for which the solution (with \code{delta=NULL}) is diagonal.
#'
#' @param Sighat The sample covariance matrix
#' @param nlam Number of lambda values to include in grid.
#' @param flmin Ratio between the smallest lambda and largest lambda in grid. (Default: 0.01)
#' Decreasing this gives less sparse solutions.
#' @param lamlist A grid of lambda values to use.  If this is non-NULL, then \code{nlam} and \code{flmin} are ignored.
#' @param w \code{(p-1)}-by-\code{(p-1)} lower-triangular matrix (above diagonal ignored).
#'          \code{w[l,]} gives the \code{l} weights for g_l.
#'          Defaults to \code{w[l,m]=sqrt(2 * l)/(l - m + 1)} for \code{m <= l}
#' @param delta Lower bound on eigenvalues.  If this is NULL (which is default), then no eigenvalue
#' constraint is included.
#' @param maxiter Number of iterations of blockwise coordinate descent to perform.
#' @param tol Only used when \code{delta} is non-NULL.  When no eigenvalue changes by more than
#' \code{tol} in BCD, convergence is assumed.
#'
#' @return Returns a sequence of convex banded estimates of the covariance matrix.
#' \describe{
#' \item{P: }{A \code{nrow(Sighat)}-by-\code{nrow(Sighat)}-by-\code{nlam} array where \code{P[, , i]} gives the \code{i}th estimate of the covariance matrix.}
#' \item{lamlist: }{Grid of lambda values used.}
#' \item{w: }{Value of w used.}
#' \item{delta: }{Value of delta used.}
#' }
#' 
#' @export
#'
#' @seealso \code{\link{hierband}} \code{\link{hierband.cv}}
#' @examples
#' set.seed(123)
#' p <- 100
#' n <- 50
#' K <- 10
#' true <- ma(p, K)
#' x <- matrix(rnorm(n*p), n, p) %*% true$A
#' Sighat <- cov(x)
#' path <- hierband.path(Sighat)
#' cv <- hierband.cv(path, x)
#' fit <- hierband(Sighat, lam=cv$lam.best)
hierband.path <- function(Sighat, nlam=20, flmin=0.01, lamlist=NULL, w=NULL, delta=NULL, maxiter=100, tol=1e-7) {
  if (!is.null(delta)) stopifnot(length(delta)==1)
  p <- nrow(Sighat)
  if (is.null(w)) {
    #w <- matrix(NA, p-1, p-1)
    #for (l in seq(p-1)) for (m in seq(l)) w[l, m] <- sqrt(2 * l) / (l - m + 1)
    w <- formw(p)
  }
  if (is.null(lamlist)) {
    lammax <- lam.max.hierband(Sighat, diag(w))
    lamlist <- lammax * exp(seq(0, log(flmin), length=nlam))
  } else {
    nlam <- length(lamlist)
  }
  P <- array(NA, c(p, p, nlam))
  for (i in seq(nlam)) {
    cat(i)
    P[, , i] <- hierband(Sighat=Sighat, lam=lamlist[i], w=w, maxiter=maxiter,
                         delta=delta, tol=tol)
  }
  cat(fill=TRUE)
  list(P=P, lamlist=lamlist, w=w, delta=delta)
}

#' Computes the smallest lambda such that P=0.
#'
#' Computes lambda_max, which is the smallest value of \code{lam} for which
#' \code{\link{hierband}} (with \code{delta=NULL}) gives a diagonal covariance matrix.
#'
#' @param Sighat empirical covariance matrix
#' @param ww the diagonal of w
lam.max.hierband <- function(Sighat, ww) {
  p <- nrow(Sighat)
  #ii <- toeplitz(seq(p, 1))
  #sighat <- rep(NA, p-1)
  #for (l in seq(p-1)) sighat[l] <- sqrt(sum(Sighat[ii==l]^2))
  sighat <- sqrt(2) * subdiagonal.l2norms(Sighat)[-1]
  max(sighat / rev(ww))
}
