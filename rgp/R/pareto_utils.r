## pareto_utils.R
##   - Utility functions for Pareto fronts 
##
## RGP - a GP system for R
## 2010-2014 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Patrick Koch, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Find the knee of a two dimensional pareto front 
##'
##' Given a matrix \code{m} of two rows and n columns, representing solutions of a
##' two-dimensional optimization problem, returns the column index of the point with
##' minimum euclidean distance to the utopia point. The utopia point is the point
##' consisting of the row minima of \code{m}. \code{NA} or \code{NaN} values of
##' \code{m} are ommited.
##'
##' @param m A matrix of two rows and n columns, representing the solutions of a
##'   two-dimensional optimization problem.
##' @param normalize Whether to normalize both objectives to the interval of
##'   [0, 1], defaults to \code{TRUE}.
##' @return The knee point index, i.e. the column index in m of the point of minimum
##'   euclidean distance to the utopia point.
##' @examples
##' m1 <- matrix(runif(200), ncol = 100)
##' plot(t(m1))
##' points(t(m1[,emoa::nds_rank(m1) == 1]), col = "red", pch = 16)
##' pKnee <- m1[, paretoFrontKneeIndex(m1)]
##' points(t(pKnee), col = "green4", pch = 16)
##'
##' @export
paretoFrontKneeIndex <- function(m, normalize = TRUE) {
  mNaOmit <- t(na.omit(t(m)))

  if (ncol(mNaOmit) == 0) {
    warning("paretoFrontKneeIndex: not enough non-NA columns in m, returning NA")
    browser() # TODO
    return (NA)
  } else if (ncol(mNaOmit) == 1) {
    return (1)
  } else {
    mNorm <- if (normalize) t(apply(mNaOmit, 1, function(row) (row - min(row)) / (max(row) - min(row)))) else mNaOmit
    pUtopia <- c(min(mNorm[1, ]), min(mNorm[2, ]))
    mDist <- as.matrix(dist(t(cbind(pUtopia, mNorm)), method = "euclidean"))[-1, 1]
    mDistMinIndex <- which.min(mDist)
    return (if (length(mDistMinIndex) == 0) 1 else mDistMinIndex) # which.min can return an empty vector
  }
}
