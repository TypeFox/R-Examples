#' Exponential p-value weights
#'
#' Computes the exponential p-value weights for multiple testing.
#' Given estimated means \code{mu} of test statistics \code{T},
#' the p-value weights are proportional to \code{exp(beta*mu)},
#' for a tilt parameter \code{beta}. In addition, the large weights are truncated
#' at a maximum value \code{UB} (upper bound), and the remaining weight is re-distributed
#' among the rest of the statistics.
#'
#' Specifically, it is assumed that \code{T} are Gaussian with mean
#' \code{mu}. One-sided tests of \code{mu>=0} against \code{mu<0}
#'  are conducted using the test statistics \code{T}. To optimize power,
#'  different levels are allocated to different tests.

#'  For more details, see the paper "Optimal Multiple Testing Under a
#'  Gaussian Prior on the Effect Sizes", by Dobriban, Fortney, Kim and Owen,
#'   \url{http://arxiv.org/abs/1504.02935}

#' @param  mu  the estimated means of the test statistics
#' @param  beta (optional) weights are proportional to \code{exp(mu*beta)}, default \code{beta=2}
#' @param  UB (optional) upper bound on the weights (default \code{UB = Inf})
#'
#' @return  The exponential weights.
#'
#' @examples
#' J <- 100
#' mu <- rnorm(J)
#' beta <- 2
#' UB <- 20
#' w <- exp_weights(mu, beta, UB)
#' @family p-value weighting
#' @seealso \code{\link{bayes_weights}} for Bayes, \code{\link{spjotvoll_weights}} for Spjotvoll weights, and \code{\link{exp_weights}} for exponential
#'   weights
#'
#' @export
#'
exp_weights <- function(mu,  beta=2,  UB = Inf) {

  #Error checking: stop if the variables are not in range
  if (UB <= 1) {
    stop("The upper bound UB must be greater than 1")
  }


  J <-  length(mu)

  #The exponential weights themselves are easy to compute
  mu_s <-  sort(mu)
  sort_index <-  order(mu)
  u <-  exp(beta * mu_s)
  c <-  mean(u)
  w <-  u / c

  #If no upper bound is set (UB==Inf), then we return the weights computed abvove. If there is a finite upper bound, then we must truncate the weights.

  if (UB == Inf) {
    #Case 1: no upper bound
    return(w)
  } else {
    #Case 2: finite upper bound
    #The previous version was written using the parametrization q = 1/UB
    #So this was not changed in the current version for simplicity
    q <- 1 / UB
    #Find large weights
    ind <- w > 1 / q

    #Find how much weight needs to be re-ditributed
    surplus <-  sum(w[ind]) - length(w[ind]) / q
    #Truncate weights (1/q = UB)
    w[ind] <-  1 / q
    k <-  length(w[ind])

    #Re-distribute the surplus among the next largest weights
    #surplus keeps track of how much extra weight we have left
    while (surplus > 0) {
      increment <-  min(1 / q - w[J - k], surplus)
      #increase the next remaining weight to the maximal possible value
      #without going above UB=1/q
      w[J - k] <-  w[J - k]  +  increment
      surplus <-  surplus  -  increment
      k <-  k + 1
    }

    #un-sort the weights
    v <-  rep(0,J)
    for (i in 1:J) {
      v[i] <-  w[sort_index == i]
    }
    w <-  v
    return(w)
  }
}
