"summary.drm" <-
  function (object, ..., digits = max(2, .Options$digits - 4), 
            correlation = TRUE) 
{
  stderr <- sqrt(diag(object$cov.scaled))
  coef <- cbind(object$coefficients, stderr, object$coefficients/stderr)
  dimnames(coef) <- list(names(object$coefficients), c("Value", 
                                                       "Std. Error", "z value"))
  object$coefficients <- coef
  if (correlation) {
    if (nrow(coef) == 1) 
      i <- 1/stderr
    else (i <- diag(1/stderr))
    object$correl <- i %*% object$cov.scaled %*% i
    dimnames(object$correl) <- list(dimnames(object$coefficients)[[1]], 
                                    dimnames(object$coefficients)[[1]])
  }
  else (object$correl <- NULL)
  object$digits <- digits
  class(object) <- "summary.drm"
  object
}
