#' @title Summarizing Gravity Model Fits
#' @description Summary method for class "gravity".
#' @param object  Object of class gravity
#' @param ...     Ignored
#' @note Summary of lme or lm gravity model, AIC and Root Mean Square Error (RMSE) of observed verses predicted 
#' @method summary gravity
#' @export
summary.gravity <- function(object, ...) {
  cat("Gravity model summary\n\n") 
  rmse <- function(y, x) { sqrt(mean((y - x)^2)) }
  print(object$formula) 
  print(summary(object$gravity))
  print( paste("AIC = ", round(object$AIC,3), sep=""))
  print( paste("RMSE = ", round(rmse(object$y, object$fit),4), sep=""))  
}
