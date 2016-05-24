## itempar generic
itempar <- function (object, ...) {
  UseMethod("itempar")
}


## methods for class 'itempar'
print.itempar <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  mdl <- attr(x, "model")
  if (mdl == "btmodel") {
    cat("Bradley-Terry model parameters (", attr(x, "model"), "):\n", sep = "")    
  } else {
    cat("Item response item parameters (", attr(x, "model"), "):\n", sep = "")
  }
  print(coef(x), digits = digits, ...)
  invisible(x)
}

coef.itempar <- function (object, ...) {
  ## remove all attributes beside names, then return named item parameters
  lbs <- names(object)
  object <- as.vector(object)
  names(object) <- lbs
  return(object)
}
    
vcov.itempar <- function (object, ...) attr(object, "vcov")

summary.itempar <- function (object, ...) {
  ## extract relevant informations
  cf <- coef(object)
  vc <- vcov(object)
  lbs <- names(object)
  
  ## compute Wald coefficient intervals for each coefficient
  citab <- cbind(cf + qnorm(p = 0.025) * sqrt(diag(vc)), cf + qnorm(p = 0.975) * sqrt(diag(vc)))
  colnames(citab) <- c("2.5 %", "97.5 %")
  rownames(citab) <- lbs

  attr(object, "summary") <- citab
  class(object) <- "summary.itempar"
  return(object)
}

print.summary.itempar <- function (x, digits = max(3, getOption("digits") - 3), 
                                   signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nItem location parameters (", attr(x, "model"), "):\n\n", sep = "")
  print(attr(x, "summary"), digits = digits, na.print = "NA", ...)
  invisible(x)
}
