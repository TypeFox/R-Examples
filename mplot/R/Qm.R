#' A measure of lack of fit
#'
#' This function calculates a lack-of-fit measure for
#' a model, \eqn{\hat{Q}_M}.  It currently simply
#' returns the negative log-likelihood value of an
#' estimated model object.
#'
#' @param object Any object from which a log-likelihood value can be extracted.
#' @param method method the model selection method to be used. Currently
#'   only \code{method = "ML"} is supported (perhaps in the future
#'   \code{method = "MVC"} will be implemented).
#' @noRd
Qm = function(object,method){
  if(method=="ML"){
    return(-as.numeric(stats::logLik(object)))
  }
}
