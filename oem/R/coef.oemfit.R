coef.oemfit <- function(object, s = NULL, ...) {
  predict(object, s = s, type = "coefficients")
}
