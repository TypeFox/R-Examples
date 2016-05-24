#' @family arc-fs
#' 
#' @template -Reweighter
#' 
#' @title Reweighter for the arc-fs algorithm.
#' @description A slightly modified implementation of the reweighter described in
#' the arc-fs algorithm.
#'
#' @details As per Leo Breiman's suggestions, a slight modification to the arc-fs
#' algorithm has been made in the case where \eqn{\epsilon}{\epsilon} -- the misclassification
#' measure -- exceed 0.5, or becomes 0. Should this happen \eqn{\beta}{\beta} is
#' set to \eqn{-\infty}{-\infty} and  a warning is produced. At this point
#' you are advised to restart the algorithm with equal probabilities or stop
#' boosting, at that iteration.
#' 
#' @return 
#' \item{beta}{scalar weights to be used by \code{\link{arcfsAggregator}}.}
#' 
#' @export
arcfsReweighter <- function(prediction, response, weights, ...) {

  d <- as.numeric(prediction != response)
  
  eps <- sum(weights * d)
  
  beta <- (1 - eps) / eps
  
  newWeights <- weights * (beta^d) / sum( weights * (beta^d) )
  
  if ( eps >= 1/2 || eps == 0 ) {
    beta <- -Inf
    newWeights <- rep.int(1/length(weights), length(weights))
    warning("beta is -inf, you'll want to remove the associated estimator.")
  }
  
  list(weights=newWeights, beta=beta)
}

class(arcfsReweighter) <- c("reweighter", class(arcfsReweighter))