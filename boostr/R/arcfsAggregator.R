#' @family arc-fs
#' @template -Aggregator
#' 
#' @title Aggregator for the arc-fs algorithm.
#' @description A (parallelized) implementation of the aggregator described in
#' the arc-fs algorithm.
#' 
#' @param beta a vector of scalar weights associated to each estimator in 
#' \code{estimators}
#' @param .parallelTally a boolean indicating if vote tallying should be
#' performed in parallel. Unless you have more than 1,000 votes / observation, 
#' you probably won't see much performance gain by parallelizing this step.
#' @param .rngSeed the RNG seed sent to
#' \code{\link{predictClassFromWeightedVote}} in the case of a tie.
#' 
#' @details By default, this function will perform its predictions in sequence
#' across the estimators in \code{estimators}. To predict in parallel, change
#' \code{.parallelPredict} to \code{TRUE}.
#' 
#' @note In accord with the arc-fs algorithm, there is the assumption that the
#' estimators in \code{estimators} are classifiers. More aptly, their output is
#' either of factor or character-type.
#' 
#' 
#' @export

arcfsAggregator <- function(estimators, beta, ..., .parallelPredict=FALSE,
                            .parallelTally=FALSE, .rngSeed=1234) {
  weights <- log(as.numeric(beta)) # cast beta as numeric, just in case.
  weightedAggregator(estimators=estimators,
                     weights=weights, 
                     .parallelPredict=.parallelPredict,
                     .parallelTally=.parallelTally,
                     .rngSeed=.rngSeed)
}

class(arcfsAggregator) <- c("aggregator", class(arcfsAggregator))