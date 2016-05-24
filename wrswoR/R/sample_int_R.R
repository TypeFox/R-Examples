#' Weighted sampling without replacement
#'
#' These functions implement weighted sampling without replacement using various
#' algorithms, i.e., they take a sample of the specified
#' \code{size} from the elements of \code{1:n} without replacement, using the
#' weights defined by \code{prob}.  The call
#' \code{sample_int_*(n, size, prob)} is equivalent
#' to \code{sample.int(n, size, replace=F, prob)}.  (The results will
#' most probably be different for the same random seed, but the
#' returned samples are distributed identically for both calls.)
#' Except for \code{sample_int_R} (which
#' has quadratic complexity as of this writing), all functions have complexity
#' \eqn{O(n \log n)}{O(n log n)} or better and
#' often run faster than R's implementation, especially when \code{n} and
#' \code{size} are large.
#'
#' @details
#'   \code{sample_int_R} is a simple wrapper for \code{\link[base]{sample.int}}.
#'
#' @inheritParams base::sample.int
#' @return An integer vector of length \code{size} with elements from
#'   \code{1:n}.
#' @seealso \code{\link[base]{sample.int}}
#' @references \url{http://stackoverflow.com/q/15113650/946850}
#' @examples
#' # Base R implementation
#' s <- sample_int_R(2000, 1000, runif(2000))
#' stopifnot(unique(s) == s)
#' p <- c(995, rep(1, 5))
#' n <- 1000
#' set.seed(42)
#' tbl <- table(replicate(sample_int_R(6, 3, p),
#'                        n = n)) / n
#' stopifnot(abs(tbl - c(1, rep(0.4, 5))) < 0.04)
#'
#' @rdname sample_int
#' @aliases sample_int_R
sample_int_R <- function(n, size, prob) {
  sample.int(n, size, replace = FALSE, prob)
}
