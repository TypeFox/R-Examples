################################################################################

predict.cv.alpha.relaxnet <- function(object,
                                      newx,
                                      alpha.val = object$which.alpha.min,
                                      type = c("link", "response",
                                               "coefficients",
                                               "nonzero", "class"),
                                      ...) {

  if(length(alpha.val) != 1)
    stop("alpha.val must be a single value")
  
  alpha.index <- which(object$alpha == alpha.val)

  if(length(alpha.index) == 0)
    stop("The specified alpha.val is not part of this cv.alpha.relaxnet object")
  
  predict(object$cv.relaxnet.results[[alpha.index]],
          newx = newx,
          type = type,
          ...)
}
  
