#' @include warehouse.R
{}



#' Return subsets of \code{AlgorithmPerformance} objects
#'
#' @param x An \code{\link{AlgorithmPerformance}} object
#' @param subset Logical expression indicating rows to keep
#' @param ... Passed to the underlying \code{\link{subset.data.frame}}
#'   call
#'
#' @return An \code{\link{AlgorithmPerformance}} object with just the
#'   selected observations
#'
#' @method subset AlgorithmPerformance
#'
#' @S3method subset AlgorithmPerformance
subset.AlgorithmPerformance <- function(x, subset,  ...) {
  e <- substitute(subset)
  r <- eval(e, x, parent.frame())
  if (!is.logical(r))
    stop("'subset' must evaluate to logical")
  r <- r & !is.na(r)

  y <- x[r, ]
  y$datasets <- y$datasets[, drop = TRUE]
  y$algorithms <- y$algorithms[, drop = TRUE]
  y$performances <- y$performances[, drop = TRUE]
  y$samples <- y$samples[, drop = TRUE]

  attr(y, "algorithm_colors") <- attr(x, "algorithm_colors")
  attr(y, "class") <- attr(x, "class")

  y
}

