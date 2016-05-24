predict.cv.sparseSVM <- function(object, X, lambda = object$lambda.min, 
                                 type=c("class","coefficients","nvars"), 
                                 exact = FALSE, ...) {
  type = match.arg(type)
  predict.sparseSVM(object$fit, X = X, lambda = lambda, 
                    type = type, exact = exact, ...)
}

coef.cv.sparseSVM <- function(object, lambda = object$lambda.min, exact = FALSE, ...) {
  coef.sparseSVM(object$fit, lambda = lambda, exact = exact, ...)
}
