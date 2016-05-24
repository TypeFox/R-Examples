print.summary.dglm <- function(x, ..., digits = NULL, quote = TRUE, prefix = "", residuals = FALSE)
{
  #  Print summary of double glm
  #  GKS  7 Jan 98
  #
  xd <- x$dispersion.summary
  x$dispersion.summary <- NULL
  if (is.null(digits))
    digits <- options()$digits
  else {
    old.digits <- options(digits = digits)
    on.exit(options(old.digits))
  }
  cat("\nCall: ")
  print(x$call)
  #
  #  Mean submodel
  #
  #   cat("\nMEAN MODEL")
  nas <- x$nas
  coef <- x$coef
  correl <- x$correl
  if (any(nas)) {
    nc <- length(nas)
    cnames <- names(nas)
    coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
    coef1[!nas,  ] <- coef
    coef <- coef1
    if (!is.null(correl)) {
      correl1 <- matrix(NA, nc, nc, dimnames = list(cnames, cnames)
      )
      correl1[!nas, !nas] <- correl
      correl <- correl1
    }
  }
  dresid <- x$deviance.resid
  n <- length(dresid)
  rdf <- x$df[2]
  if (residuals) {
    if (rdf > 5) {
      cat("Deviance Residuals:\n")
      rq <- stats::quantile(as.vector(dresid))
      names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
      print(rq, digits = digits)
    }
    else if (rdf > 0) {
      cat("Deviance Residuals:\n")
      print(dresid, digits = digits)
    }
  }
  if (any(nas))
    cat("\nMean Coefficients: (", sum(nas), 
        " not defined because of singularities)\n", sep = "")
  else cat("\nMean Coefficients:\n")
  print(coef, digits = digits)
  
  cat(paste("(Dispersion Parameters for", x$family$family, 
            "family estimated as below", ")\n"))
  int <- attr(x$terms, "intercept")
  if (is.null(int))
    int <- 1
  cat("\n    Scaled Null Deviance:", format(round(x$null.deviance, digits)), "on", n - 
        int, "degrees of freedom\n")
  cat("Scaled Residual Deviance:", format(round(x$deviance, digits)), "on", round(
    rdf, digits), "degrees of freedom\n")
  #       cat("\nNumber of Fisher Scoring Iterations:", format(trunc(x$iter)), "\n")
  if (!is.null(correl)) {
    p <- dim(correl)[2]
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1,  -p, drop = FALSE], quote = FALSE, digits = digits)
    }
  }
  #
  #  Dispersion submodel
  #
  nas <- xd$nas
  coef <- xd$coef
  correl <- xd$correl
  if (any(nas)) {
    nc <- length(nas)
    cnames <- names(nas)
    coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
    coef1[!nas,  ] <- coef
    coef <- coef1
    if (!is.null(correl)) {
      correl1 <- matrix(NA, nc, nc, dimnames = list(cnames, cnames)
      )
      correl1[!nas, !nas] <- correl
      correl <- correl1
    }
  }
  dresid <- xd$deviance.resid
  n <- length(dresid)
  rdf <- xd$df[2]
  if (residuals) {
    if (rdf > 5) {
      cat("Deviance Residuals:\n")
      rq <- stats::quantile(as.vector(dresid))
      names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
      print(rq, digits = digits)
    }
    else if (rdf > 0) {
      cat("Deviance Residuals:\n")
      print(dresid, digits = digits)
    }
  }
  if (any(nas))
    cat("\nDispersion Coefficients: (", sum(nas), 
        " not defined because of singularities)\n", sep = "")
  else cat("\nDispersion Coefficients:\n")
  print(coef, digits = digits)
  cat(paste("(Dispersion parameter for", xd$family$family, 
            "family taken to be", format(round(xd$dispersion, digits)), ")\n"))
  int <- attr(xd$terms, "intercept")
  if (is.null(int))
    int <- 1
  cat("\n    Scaled Null Deviance:", format(round(xd$null.deviance, digits)), "on", n - 
        int, "degrees of freedom\n")
  cat("Scaled Residual Deviance:", format(round(xd$deviance, digits)), "on", round(
    rdf, digits), "degrees of freedom\n")
  #   cat("\nNumber of Fisher Scoring Iterations:", format(trunc(xd$iter)), "\n")
  if (!is.null(correl)) {
    p <- dim(correl)[2]
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1,  -p, drop = FALSE], quote = FALSE, digits = digits)
    }
  }
  #
  #  Overall iteration
  #
  cat("\nMinus Twice the Log-Likelihood:", format(round(x$m2loglik, digits)), "\n")
  cat("Number of Alternating Iterations:", format(trunc(x$outer.iter)), "\n")
  invisible(NULL)
}
