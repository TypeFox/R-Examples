#' @family Wrapper Generators
#' 
#' @title Create a boostr compatible wrapper for an aggregator.
#' @description Use provided metadata on a given aggregator to create a boostr
#' compatible wrapper. See section below for details on aggregators.
#' 
#' @param aggregator a function which satisfies the abstract definition of an 
#' aggregator.
#' @param .inputEnsemble a string indicating the name of the argument that
#' \code{aggregator} uses for the ensemble of estimators created during the
#' Boost algorithm.
#' @param .verbose a logical flag indicating whether warnings should be output
#' or not.
#' 
#' @template aggregators
#' 
#' @return A function with is also an '\code{aggregator}' object. The function's
#' signature and output are now compatible with the boostr framework. In
#' particular, the signature of the wrapper is 
#' \item{estimators}{the list of estimators to be sent to \code{aggregator}.}
#' \item{...}{any additional arguments accepted/required by \code{aggregator}.}
#' The output of this aggregator is an estimator with signature
#' \item{newdata}{the data the \code{aggregator}'s output would use for prediction.}
#' 
#' @examples
#' \dontrun{
#' testAggregator <- function(ensemble) {
#'  weights <- runif(min=0, max=1, n=length(ensemble))
#'  function(x) {
#'    preds <- foreach(estimator = iter(ensemble),
#'                   .combine = rbind) %do% {
#'                     matrix(as.character(estimator(x)), nrow=1)
#'                   }
#'  
#'    as.factor(predictClassFromWeightedVote(preds, weights))
#'  }
#' }
#'
#' wrappedAggregator <- wrapAggregator(testAggregator,
#'                                     .inputEnsemble="ensemble")
#'}                                     


wrapAggregator <- function(aggregator, .inputEnsemble="estimators", .verbose=FALSE) {

  # extend aggregator signature.
  aggregator <- addDots(aggregator, .verbose)
  
  f <- function(estimators, ...) {
    aggregatorArgs <- c(list(estimators=estimators), list(...))
    names(aggregatorArgs)[1] <- .inputEnsemble
    
    estimator <- do.call(aggregator, aggregatorArgs)
    
    if (length(names(formals(estimator))) > 1 ) {
      stop('resulting estimator from aggregator should only take one argument')
    }
    
    # estimator should be a univariate function
    function(newdata) {
      estimator(newdata)
    }
  }
  class(f) <- c("aggregator", class(f))
  f
}

