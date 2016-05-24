summary.aodql <- function(object, ...) {
  Z <- switch(object$method, chisq = "X2", dev = "D")
	cat("\nMu coefficients:\n")	
	print(summary(object$fm))
  cat("Method: ", object$method, "\n", sep = "")
  cat("Phi coefficient: ", object$phi, "\n", sep = "")
  cat(Z, " = ",  object$Z, " (nbiter = ",
    object$nbiter[1], "/", object$nbiter[2], ")", sep = "")
  }
