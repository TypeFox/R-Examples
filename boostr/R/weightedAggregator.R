#' @param weights a vector of scalar weights associated to each estimator in 
#' \code{estimators}
#' 
#' @details If performing regression and your estimators produce \code{NA}'s, you
#' can have \code{\link{weighted.mean}} remove the \code{NA}'s by passing
#' \code{na.rm=TRUE} to \code{weightedAggregator}'s function call.
#' 
#' @template -Aggregator
#' @template stockAggregators
#' 
#' @rdname stockAggregators
#' @export

weightedAggregator <- function(estimators, weights, ..., .parallelPredict=FALSE,
                            .parallelTally=FALSE, .rngSeed=1234) {
  
  weights <- as.numeric(weights)
  
  function(newdata) {
    preds <- makePredictions(estimators, newdata, .parallelPredict)
    
    if (typeof(preds) == "character") {
      out <- do.call(predictClassFromWeightedVote, 
                     list(preds=preds, weights=weights, .parallel=.parallelTally,
                          .rngSeed=.rngSeed))
      
      as.factor(out)
    } else {
      do.call(predictResponseFromWeightedAverage,
              list(preds=preds, weights=weights,
                   .parralel=.parallelTally, list(...)))
    }
  }
}

class(weightedAggregator) <- c("aggregator", class(weightedAggregator))