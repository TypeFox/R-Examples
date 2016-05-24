#' @title
#' Standard (vanilla) bagging procedure.
#' @description
#' Build an estimator from a simple resampling of data.
#' 
#' @param prediction a vector of predictions.
#' @param response a vector whose \eqn{i^{th}}{ith} component is the true
#' response for the \eqn{i^{th}}{ith} component of \code{prediction}.
#' @param ... implemented to allow reweighter to accept its output as its input.
#' 
#' @note a "bagger" is just a reweighter who returns uniform weights regardless of
#' the input.
#' 
#' @export
#' 
#' @return a list with component '\code{weights}': 
#' a normalized vector of 1's with length equal to that of \code{response}.
#' 
#' @family reweighters
#'

vanillaBagger <- function(prediction, response, ...){
  list(weights=rep.int(1, length(prediction)))
}

class(vanillaBagger) <- c("reweighter", class(vanillaBagger))