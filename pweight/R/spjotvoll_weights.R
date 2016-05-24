#' Spjotvoll p-value weights
#'
#' Computes the Spjotvoll p-value weights for multiple testing.
#' Given estimated means \code{mu} of test statistics \code{T},
#' the weighting scheme optimizes the expected number of discoveries
#' using \code{T} as test statistics in multiple testing
#' at some specific level \code{q}.
#'
#' Specifically, it is assumed that \code{T} are Gaussian with mean
#' \code{mu}. One-sided tests of \code{mu>=0} against \code{mu<0}
#'  are conducted using the test statistics \code{T}. To optimize power,
#'  different levels are allocated to different tests.

#'  For more details, see the paper "Optimal Multiple Testing Under a
#'  Gaussian Prior on the Effect Sizes", by Dobriban, Fortney, Kim and Owen,
#'   \url{http://arxiv.org/abs/1504.02935}
#'
#'@param mu a negative vector, the estimated means of test statistics
#'@param q level at which tests will be performed
#
#'@return The optimal Spjotvoll weights.
#'
#' @examples
#' J <- 100
#' mu <- -abs(rnorm(J))
#' q <- 0.05 / J
#' w <- spjotvoll_weights(mu, q)
#'
#' @family p-value weighting
#'
#'
#'@export
#'
spjotvoll_weights <- function(mu, q) {

  #Error checking: stop if the variables are not in range
  if (any(mu >= 0)) {
    stop("Negative means mu required")
  }

  if ((q <= 0) |
      (q >= 1)) {
    stop("Level at which tests will be performed must be in (0,1)")
  }

  J <- length(mu)

  #define the functions used by Newton's method
  #the goal is to find the zero of the function f below
  f <- function(c)
    1 / J * sum(pnorm(c / mu + mu / 2)) - q;
  df <- function(c)
    1 / J * sum(dnorm(c / mu + mu / 2) / mu);

  # set parameters of Newton iteration:
  # starting point, tolerance, max number of iterations
  x0 <- 0;
  tol <- 1e-3 / J;
  nmax <- 100;
  #these vectors will store the trajectory of the iteration,and the step sizes
  x <- rep(0, nmax)
  ex <- rep(0, nmax)

  #Newton's method iterations
  #set first step
  x[1] <- x0 - f(x0) / df(x0);
  ex[1] <- abs(x[1] - x0);
  k <- 2;

  #iterate while step size is sufficiently large
  while ((ex[k - 1] >= tol) && (k <= nmax)) {
    x[k] <- x[k - 1] - (f(x[k - 1]) / df(x[k - 1]));
    ex[k] <- abs(x[k] - x[k - 1]);
    k <- k + 1;
  }

  #the last step of the trajectory
  c <- x[k - 1]; #could return this

  #compute weights
  w <- pnorm(c / mu + mu / 2) / q;

  return(w)
}
