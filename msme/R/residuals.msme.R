residuals.msme <- function(object, type = c("deviance","standard"), ...) {
  type <- match.arg(type)
  if (type == "standard") {
    object$residuals / sqrt(1 - diag(hatvalues(object)))
  } else {
    object$residuals
  }
}


