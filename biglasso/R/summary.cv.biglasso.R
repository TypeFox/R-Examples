summary.cv.biglasso <- function(object, ...) {
  # inherits cv.ncvreg
  class(object) <- 'cv.ncvreg'
  summary(object = object, ...)
}

print.summary.cv.biglasso <- function(x, digits, ...) {
  # inherits summary.cv.ncvreg
  class(x) <- 'summary.cv.ncvreg'
  print(x = x, ...)
}