#'
#' @rdname predictClassViaVoting
#' @aliases predictClassFromWeightedVote
#' @title Predict a class using (un)weighted voting.
#' @description Process a matrix of class predictions and form a column-wise
#' estimate based on weighted voting.
#'
#'  
#' @param preds is (character) matrix of predicted classes
#' @param weights is a vector of length equal to \code{nrow(preds)}
#' @param .parallel is a boolean flag determining whether to work
#' across columns of \code{preds} in parallel -- need to register a
#' parallel backend (e.g. \code{doParallel}, \code{doRedis}) for this to
#' actually work.
#' @param .rngSeed the value of the RNG seed to be used in the case that 
#' ties are to be randomly broken. 
#' 
#' @details Gives the vote from row(i) in \code{preds} weight equal to
#' \code{weights[i]}. Ties are broken randomly, but before so, the seed is set
#' to \code{.rngSeed}. 
#' 
#' @export
#' 
#' @return a character vector of length equal to \code{ncol(preds)} containing
#' the class estimates per column of \code{preds}. 

predictClassFromWeightedVote <- function(preds, weights, .parallel=FALSE,
                                         .rngSeed=1234) {
  determineWinner <- function(votes) {
    # entries of votes represent individual levels
    winners <- which(votes == max(votes))
    if (length(winners) > 1) {
      set.seed(.rngSeed)
      sample(names(votes)[winners], size=1)
    } else {
      names(votes)[winners]
    }
  }

  `%op%` <- if (getDoParRegistered() && .parallel) `%dopar%` else `%do%`
  
  # may want to consider making this step parallel
  # my gut says that its too fast to benefit, though.
  
  prediction <- NULL # instantiate local variable
  
  foreach(prediction = iter(preds, by='column'), .combine=c) %op% {
    # tally over the unique responses in prediction column
    votes <- lapply(unique(prediction), function(level) {
      # weight the predictions
      out <- as.vector(weights %*% as.numeric(prediction == level))
      names(out) <- level
      return(out)
    })
    determineWinner(unlist(votes))
  }
}

#' @aliases predictClassFromVote
#' @export
#' @rdname predictClassViaVoting
#' 
predictClassFromVote <- function(preds, .parallel=FALSE, .rngSeed=1234) {
  if (!is.null(dim(preds))) {
    weights <- rep.int(1, nrow(preds))
  } else {
    weights <- rep.int(1, length(preds))
  }
  predictClassFromWeightedVote(preds, weights, .parallel=.parallel,
                               .rngSeed=.rngSeed)
}