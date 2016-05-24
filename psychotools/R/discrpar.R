## discrpar generic
discrpar <- function (object, ...) {
  UseMethod("discrpar")
}


## methods for class 'discrpar'
print.discrpar <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Item response discrimination parameters (", attr(x, "model"), "):\n", sep = "")
  print(coef(x), digits = digits, ...)
  invisible(x)
}

coef.discrpar <- function (object, ...) {
  ## remove all attributes beside names, then return named item parameters
  lbs <- names(object)
  object <- as.vector(object)
  names(object) <- lbs
  return(object)
}
    
vcov.discrpar <- function (object, ...) attr(object, "vcov")
