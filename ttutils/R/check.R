check <- function(object, ...) {
  UseMethod("check")
}

check.default <- function(object, ...) {
  warning("No 'check' function defined for class '", class(object), "'!")
  return(FALSE)
}
