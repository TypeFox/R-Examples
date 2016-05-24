#' @aliases arcx4Aggregator
#' @details
#' \code{arcx4Aggregator} is just \code{vanillaAggregator} by another name.
#' @rdname stockAggregators
#' @family arcx4
#' @export
arcx4Aggregator <- function(estimators, ..., .parallelPredict=FALSE,
                            .parallelTally=FALSE, .rngSeed=1234) {
  vanillaAggregator(estimators=estimators, .parallelPredict=.parallelPredict,
                    .parallelTally=.parallelTally, .rngSeed=.rngSeed)
}

class(arcx4Aggregator) <- c("aggregator", class(arcx4Aggregator))
