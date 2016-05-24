#' @title Naive multi-start version of \code{\link{optimize}} for global optimization.
#'
#' @description
#' The univariate \code{\link{optimize}} can stop at arbitrarily bad points when
#' \code{f} is not unimodal. This functions mitigates this effect in a very naive way:
#' \code{interval} is subdivided into \code{nsub} equally sized subintervals,
#' \code{\link{optimize}} is run on all of them (and on the original big interval) and
#' the best obtained point is returned.
#'
#' @param f See \code{\link{optimize}}.
#' @param interval See \code{\link{optimize}}.
#' @param ... See \code{\link{optimize}}.
#' @param lower See \code{\link{optimize}}.
#' @param upper See \code{\link{optimize}}.
#' @param maximum See \code{\link{optimize}}.
#' @param tol See \code{\link{optimize}}.
#' @param nsub [\code{integer(1)}]\cr
#'   Number of subintervals. A value of 1 implies normal \code{\link{optimize}} behavior.
#'   Default is 50L.
#' @return See \code{\link{optimize}}.
#' @export
optimizeSubInts = function(f, interval, ..., lower = min(interval), upper = max(interval),
  maximum = FALSE, tol = .Machine$double.eps^0.25, nsub = 50L) {
  nsub = asCount(nsub, positive = TRUE)

  interval = c(lower, upper)
  # run on normal interval
  best = optimize(f = f, interval = interval, maximum = maximum, tol = tol)

  # run on smaller partitions
  if (nsub > 1L) {
    mult = ifelse(maximum, -1, 1)
    grid = seq(lower, upper, length.out = nsub - 1L)
    for (j in seq_len(length(grid)-1L)) {
      res = optimize(f = f, interval = c(grid[j], grid[j+1L]), maximum = maximum, tol = tol)
      if (mult * res$objective < mult * best$objective)
        best = res
    }
  }
  return(best)
}
