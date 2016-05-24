`Cp.lm` <-
function(object, a = 1, b = 0, S2, ...){
  a * (crossprod(object$residuals, object$residuals)/S2 - length(object$residuals) +
  2*length(object$coefficients)) + b * (abs(crossprod(object$residuals, object$residuals)/S2 -
  length(object$residuals) + length(object$coefficients)))
  }

