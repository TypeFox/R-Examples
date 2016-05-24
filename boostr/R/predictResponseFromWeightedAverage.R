#' @title Predict a numeric response using (un)weighted averaging.
#' @description Process a matrix of predicted responses and form a column-wise
#' estimate based on (un)weighted averaging.
#'
#'  
#' @param preds is matrix of predicted classes
#' @param weights is a vector of length equal to \code{nrow(preds)}
#' @param .parallel is a boolean flag determining whether to work
#' across columns of \code{preds} in parallel -- need to register a
#' parallel backend (e.g. \code{doParallel}, \code{doRedis}) for this to
#' actually work.
#' @param ... additional arguments to pass to \code{\link{weighted.mean}}.
#' 
#' @details Gives the prediction from row(i) in \code{preds} weight equal to
#' \code{weights[i]}. Note that \code{NA}'s are not removed. To have
#' \code{\link{weighted.mean}} remove the \code{NA}'s pass \code{na.rm=TRUE} to
#' the function call.
#' 
#' @export
#' 
#' @return a vector of length equal to \code{ncol(preds)} containing the
#' estimated response for each column of \code{preds}. 
predictResponseFromWeightedAverage <- function(preds, weights, .parallel, ...) {
  `%op%` <- if (getDoParRegistered() && .parallel) `%dopar%` else `%do%`
  
  predsVec <- NULL # instantiate local variable
  
  foreach(predsVec = iter(preds, by="column"), .combine=c) %op% {
    weighted.mean(x=predsVec, w=weights, ...)
  }
}