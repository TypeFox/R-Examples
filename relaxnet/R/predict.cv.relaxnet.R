################################################################################

predict.cv.relaxnet <- function(object,
                                newx,
                                which.model = object$which.model.min,
                                s = object$overall.lambda.min,
                                type = c("link", "response", "coefficients",
                                         "nonzero", "class"),
                                exact = FALSE,
                                ...) {

  type = match.arg(type)

  if(!identical(exact, FALSE))
    stop("In this version, only exact = FALSE has been implemented")
  
  predict(object$relaxnet.fit,
          newx,
          which.model,
          s,
          type,
          exact,
          ...)
  
}
