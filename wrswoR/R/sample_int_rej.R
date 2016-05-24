#' @rdname sample_int
#'
#' @details \code{sample_int_rej} uses repeated weighted sampling with
#'   replacement and a variant of rejection sampling. It is implemented purely
#'   in R.
#'   This function simulates weighted sampling without replacement using
#'   somewhat more draws \emph{with} replacement, and then discarding
#'   duplicate values (rejection sampling).  If too few items are
#'   sampled, the routine calls itself recursively on a (hopefully) much
#'   smaller problem.  See also
#'   \url{http://stats.stackexchange.com/q/20590/6432}.
#' @author Kirill Müller (for \code{*_int_rej})
#' @examples
#' ## Rejection sampling
#' s <- sample_int_rej(20000, 10000, runif(20000))
#' stopifnot(unique(s) == s)
#' p <- c(995, rep(1, 5))
#' n <- 1000
#' set.seed(42)
#' tbl <- table(replicate(sample_int_rej(6, 3, p),
#'                        n = n)) / n
#' stopifnot(abs(tbl - c(1, rep(0.4, 5))) < 0.04)
#'
sample_int_rej <- function(n, size, prob) {
  .check_args(n, size, prob)
  .sample_int_rej(n, size, prob, 2, 1)
}

# Euler-Mascheroni constant
.EM = 0.57721566490153286060651209008240243104215933593992

# Computes the harmonic series. Exact for the first
# .harmonic.series.max values (through table lookup), otherwise using
# the approximation ln(a) + \gamma + 1 / (2a). Source:
# http://en.wikipedia.org/wiki/Harmonic_number
.harmonic <- function(a) {
  stopifnot(a >= 0)
  if (a < length(.harmonic.series)) {
    .harmonic.series[a + 1]
  } else {
    log(a) + .EM + .5 / a
  }
}

#' @importFrom logging logdebug
# Workhorse
.sample_int_rej <- function(
  n, size, prob, MAX_OVERSHOOT, BIAS) {

  logdebug('.sample_int_rej: parameters: %s, %s, %s', n, size, length(prob))

  # How many draws *with replacement* are required on average, assuming
  # *uniform* weights? (With non-uniform weights, this number can only
  # increase.) The result is a general case of the coupon collector
  # problem, see http://math.stackexchange.com/q/247569/16420 for an
  # analysis. BIAS can be supplied to correct the estimate by a factor,
  # at most n * MAX_OVERSHOOT samples will be drawn.  Both are tuning
  # parameters, ideal values are still to be found through simulation.
  wr.size <- ceiling(n * min(BIAS * (.harmonic(n) - .harmonic(n - size)),
                             MAX_OVERSHOOT))
  logdebug('.sample_int_rej: wr.size=%s', wr.size)

  # Do the sampling with replacement...
  wr.sample <- sample.int(n, size=wr.size, replace=T, prob)
  # ...but keep only unique values.
  wr.sample <- unique(wr.sample)
  wr.sample.len <- length(wr.sample)
  logdebug('.sample_int_rej: wr.sample.len=%s', wr.sample.len)

  # How much still left to do?
  rem.size <- size - wr.sample.len
  # Done? Great!
  if (rem.size <= 0)
    return (head(wr.sample, size))

  # Not yet: Find out which indexes haven't been sampled yet.  Recall
  # that negative indexes in a vector subscription mean "all but
  # the selected".
  rem.indexes <- (1:n)[-wr.sample]
  rem.n <- length(rem.indexes)
  stopifnot(rem.n == n - wr.sample.len)

  # Recursive call to sample without replacement from the remaining
  # weights
  rem.sample <- .sample_int_rej(rem.n, rem.size,
                                prob[rem.indexes],
                                MAX_OVERSHOOT, BIAS)

  # Combine the results, substitute the indexes from 1:rem.n obtained
  # from the recursive call using the rem.indexes map
  c(wr.sample, rem.indexes[rem.sample])
}
