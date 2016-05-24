gof.default <- function(object)
  stop("No applicable method for objects of class ", class(object))

gof <- function(object){
  D <- deviance(object)
  X2  <- sum(residuals(object, type = "pearson")^2)
  structure(list(D = D, X2 = X2, df = df.residual(object)), class = "gof")
  }
