#' @rdname stockAggregators
#' 
#' @export
vanillaAggregator <- function(estimators, ..., .parallelPredict=FALSE,
                              .parallelTally=FALSE, .rngSeed=1234) {
  
  wts <- rep.int(1, length(estimators))

  weightedAggregator(estimators=estimators, weights=wts,
                     .parallelPredict=.parallelPredict,
                    .parallelTally=.parallelTally, .rngSeed=.rngSeed) 
}

class(vanillaAggregator) <- c("aggregator", class(vanillaAggregator))