#' @title Methods for Cdfqr Objects
#' @aliases predict.cdfqr fitted.cdfqr
#' @description Methods for obtaining the fitted/predicted values for a fitted cdfqr object.
#' @param object A cdfqr model fit object
#' @param type A character that indicates whether the full model prediction/fitted values are needed, or values for the `mu` and `sigma` submodel only. 
#' @param ... currently ignored
#' @examples
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2','t2', data = JurorData)
#' 
#' plot(predict(fit))
#' 
#' @method predict cdfqr 
#' @export 
#' 
predict.cdfqr <- function(object, type = c("full","mu","sigma"), ...) {
  type <- match.arg(type)
  
  fitted <- switch(type, full = {
    object$fitted$full
  }, mu = {
    object$fitted$full_mu
  }, sigma = {
    object$fitted$full_sigma
  })
  
  return(fitted)
}

#' @method fitted cdfqr 
#' @export 
#' @rdname predict.cdfqr
fitted.cdfqr <- function(object, type = c("full","mu","sigma"), ...) {
  type <- match.arg(type)
  
  fitted <- switch(type, full = {
    object$fitted$full
  }, mu = {
    object$fitted$mu
  }, sigma = {
    object$fitted$sigma
  })
  
  return(fitted)
}