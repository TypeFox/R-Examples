#' @family adaboost
#' @template -Aggregator
#' 
#' @title Aggregator for the Adaboost.M1 algorithm
#' @description Implements a (parallelized) version of the aggregator described
#' in the Adaboost.M1 algorithm.
#' 
#' @param alpha a vector (or list) of length equal to the length of 
#' \code{estimators}. Each entry of \code{alpha} acts as a prediction weight for
#' the corresponding estimator.
#' 
#' @export
adaboostAggregator <- function(estimators, alpha, ..., .parallelPredict=FALSE) {
  function(newdata) {
    # makeClassPrediction returns characters, not numbers
    preds <- makePredictions(estimators, newdata, .parallelPredict)

    # recast preds as a numeric matrix
    predDims <- dim(preds)
    preds <- matrix(as.numeric(preds), nrow=predDims[1], ncol=predDims[2])
    
    # aggregate predictions and then re-cast them as a factor variable.
    factor(sign( as.numeric(alpha) %*% preds ))
  }
}

class(adaboostAggregator) <- c("aggregator", class(adaboostAggregator))