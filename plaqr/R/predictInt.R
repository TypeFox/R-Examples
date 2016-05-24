
predictInt <- function(fit, level=.95, newdata=NULL, ...)
{
  x <- fit
  taulwr <- (1-level)/2
  tauupr <- .5+level/2

  # If newdata is NULL, use current values for prediction
  if(is.null(newdata)){

    # Median
    if(fit$tau==.5){
      median <- fit$fitted.values
    } else {
      x$call$tau <- .5
      median <- eval.parent(x$call)$fitted.values
    }

    # Lower quantile
    if(fit$tau==taulwr){
      lwr <- fit$fitted.values
    } else {
      x$call$tau <- taulwr
      lwr <- eval.parent(x$call)$fitted.values
    }

    # Upper quantile
    if(fit$tau==tauupr){
      upr <- fit$fitted.values
    } else {
      x$call$tau <- tauupr
      upr <- eval.parent(x$call)$fitted.values
    }

  } else {

    # Median
    if(fit$tau==.5){
      median <- predict(fit, newdata)
    } else {
      x$call$tau <- .5
      median <- predict(eval.parent(x$call), newdata)
    }

    # Lower quantile
    if(fit$tau==taulwr){
      lwr <- predict(fit, newdata)
    } else {
      x$call$tau <- taulwr
      lwr <- predict(eval.parent(x$call), newdata)
    }

    # Upper quantile
    if(fit$tau==tauupr){
      upr <- predict(fit, newdata)
    } else {
      x$call$tau <- tauupr
      upr <- predict(eval.parent(x$call), newdata)
    }

  }
  mat <- cbind(median,lwr,upr)
  return(mat)
}