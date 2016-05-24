#' Predict from a cosinor model
#'
#' Given a time variable and optional covariates, generate predicted values from
#' a cosinor fit. Default prediction is the mean value, optionally can predict
#' at a given month
#'
#' @param object An object of class \code{cosinor.lm}
#' @param newdata Optional new data
#' @param ... other arguments
#'
#'
#' @examples
#'
#' fit <- cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
#' predict(fit)
#'
#' @export
#'

predict.cosinor.lm <- function(object, newdata, ...){

  if(missing(newdata) || is.null(newdata)) {

  Y <- object$fit$model[,paste(attr(object$Terms, "variables")[1 + attr(object$Terms, "response")])]
  Y.hat <- fitted(object$fit)

} else {

  Y <- newdata[, paste(attr(object$Terms, "variables")[1 + attr(object$Terms, "response")])]
  Y.hat <- predict(object$fit, newdata = newdata)

}

mu.hat <- object$fit$coefficients[1]
return(Y - Y.hat + mu.hat)

}

