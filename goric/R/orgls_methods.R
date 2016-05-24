summary.orgls <- function(object, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("logLik: ", format(object$logLik, digits = digits),"\n\n")
  if (length(coef(object))) {
    cat("Coefficients:\n")
    print.default(format(coef(object), digits = digits), print.gap = 2, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  if (length(coef(object))) {
    cat("Unconstrained solution:\n")
    print.default(format(object$unccoefficients, digits = digits), print.gap = 2, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  cat("\nActive constraints:",object$iact,"\n")
  if (length(object$modelStruct) > 0) {
    print(summary(object$modelStruct))
  }
  invisible(object)
}

print.orgls <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  if (length(x$modelStruct) > 0) {
    print(summary(x$modelStruct))
  }
  invisible(x)
}


model.matrix.orgls <- function(object, ...){
  object$XW
}


predict.orgls <- function(object, newdata=NULL, ...){
  print("Not yet implemented...")
}

