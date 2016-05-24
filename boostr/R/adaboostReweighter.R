#' @family adaboost
#' @template -Reweighter
#' 
#' @title Reweighter function for the Adaboost.M1 algorithm
#' @description Implements a slightly modified version of the reweighter described
#' in the Adaboost.M1 algorithm.
#' 
#' @details The modification of the reweighter comes in to play when 
#' \eqn{$\epsilon = 0$}{\epsilon = 0}. This is when the esimator correctly
#' classifies every observation in the learning set. Consequently, we're
#' supposed to define 
#' \eqn{
#' \alpha = \displaystyle log\left(\frac{1-\epsilon}{\epsilon}\right)}{
#' \alpha = log(1-\epsilon) - log(\epsilon)
#' }
#' However, this is \eqn{+\infty}{+\infty}, which is not a number R is used
#' to working with. To work around this, we have to create a conditional statement
#' that sets \code{alpha <- log(.Machine$double.xmax)} and let the algorithm
#' proceed as originally described. The effect of this modification is the following:
#' \enumerate{
#'  \item{the update that's supposed to be made to \code{weights}, which is a 
#'  function of \code{alpha}, effectively keeps \code{weights} as they were
#'   before.}
#'  \item{if you pair this reweighter with \code{\link{adaboostAggregator}} then
#'  the estimator associated to this very large \code{alpha} now has tremendous
#'  weight inside the weighted sum in the aggregator. This isn't, necessarily,
#'  a bad thing -- the estimator classified every observation in \code{data}
#'  correctly.}
#' }
#' 
#' @return
#' \item{alpha}{performance measure of \code{estimator} to be used by
#' \code{\link{adaboostAggregator}}.}
#'
#' @export
adaboostReweighter <- function(prediction, response, weights, ...) {
  d <- as.numeric(prediction != response)
  
  err <- sum( weights * d ) / sum(weights)
  
  alpha <- if (err != 0) log((1-err) / err) else log(.Machine$double.xmax)
  
  newWeights <- weights * exp(alpha * d)
  
  list(weights=newWeights, alpha=alpha)
}

class(adaboostReweighter) <- c("reweighter", class(adaboostReweighter))