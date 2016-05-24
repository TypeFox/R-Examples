summary.orlm <- function(object, digits = max(3, getOption("digits") - 3), ...){
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
  invisible(object)
}

print.orlm <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(x)
}


model.matrix.orlm <- function(object, ...){
  object$X
}


predict.orlm <- function(object, newdata=NULL, ...){
  if (is.null(newdata)){
    X <- object$X
  } else {
    form <- object$call$formula
    form[[2]] <- NULL
    X <- model.matrix(as.formula(form), data=newdata)
  }
  X %*% object$coefficients
}

