## personpar generic
personpar <- function (object, ...) {
  UseMethod("personpar")
}


## methods for class 'personpar'
print.personpar <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Item response person parameters (", attr(x, "model"), "):\n", sep = "")
  print(coef(x), digits = digits, ...)
  invisible(x)
}

coef.personpar <- function (object, ...) {
  ## remove all attributes beside names, then return named item parameters
  lbs <- names(object)
  object <- as.vector(object)
  names(object) <- lbs
  return(object)
}

vcov.personpar <- function (object, ...) {
  ## extract covariance matrix, if not existent, create NA matrix
  vc <- attr(object, "vcov")
  if (is.null(vc)) matrix(NA, ncol = length(object), nrow = length(object)) else vc
}

summary.personpar <- function (object, ...) {
  ## extract relevant informations
  cf <- coef(object)
  vc <- vcov(object)
  lbs <- names(object)
  
  ## compute Wald coefficient intervals for each coefficient
  citab <- cbind(cf + qnorm(p = 0.025) * sqrt(diag(vc)), cf + qnorm(p = 0.975) * sqrt(diag(vc)))
  colnames(citab) <- c("2.5 %", "97.5 %")
  rownames(citab) <- lbs

  attr(object, "summary") <- citab
  class(object) <- "summary.personpar"
  return(object)
}

print.summary.personpar <- function (x, digits = max(3, getOption("digits") - 3), 
                                     signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\n Person parameters (", attr(x, "model"), "):\n\n", sep = "")
  print(attr(x, "summary"), digits = digits, na.print = "NA", ...)
  invisible(x)
}
