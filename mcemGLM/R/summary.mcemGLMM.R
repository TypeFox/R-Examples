summary.mcemGLMM <- function(object, ...) {
  # Fixed effects
  coef0 <- tail(object$mcemEST, n = 1)[1:ncol(object$x)]
  names(coef0) <- colnames(object$mcemEST)[1:ncol(object$x)]
  
  ran.eff0 <- colMeans(object$randeff)
  
  # Variance estimates
  var.est0 <-tail(object$mcemEST, n = 1)[-(1:ncol(object$x))]
  names(var.est0) <- colnames(object$mcemEST)[-(1:ncol(object$x))]
  
  # Covariance matrix and standard errors
  cmat0 <- solve(object$iMatrix)
  std.err0 <- sqrt(diag(cmat0))
  std.err1 <- std.err0[-(1:ncol(object$x))]
  std.err0 <- std.err0[1:ncol(object$x)]
  
  # z values
  zval0 <- coef0/std.err0[1:ncol(object$x)]
  zval1 <- var.est0/std.err1
  
  # p values
  pval0 <- 2 * pnorm(-abs(zval0))
  pval1 <- pnorm(-abs(zval1))
  
  resultsFixed <- matrix(0, length(coef0), 4)
  resultsFixed[, 1] <- coef0
  resultsFixed[, 2] <- std.err0[1:ncol(object$x)]
  resultsFixed[, 3] <- zval0
  resultsFixed[, 4] <- round(pval0, 8)
  rownames(resultsFixed) <- names(coef0)
  colnames(resultsFixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  resultsVar <- matrix(0, length(var.est0), 4)
  resultsVar[, 1] <- var.est0
  resultsVar[, 2] <- std.err1
  resultsVar[, 3] <- zval1
  resultsVar[, 4] <- pval1
  rownames(resultsVar) <- names(var.est0)
  colnames(resultsVar)  <- c("Estimate", "Std. Error", "z value", "Pr(>z)")
  
  cat("Call:\n  ")
  print(object$call)
  
  cat("\n   Two sided Wald tests for fixed effects coefficients:\n\n")
  print(resultsFixed)
  
  if (object$call$family %in% c("bernoulli", "poisson")) {
    cat("\n\n   One sided Wald tests for variance components:\n\n")
    print(resultsVar)
  } else {
    resultsAlpha <- matrix(resultsVar[1, 1:2], 1, 2)
    # resultsTheta <- matrix(0, 1, 2)
    # resultsTheta[1, 1] <- 1 + 1/resultsAlpha[1, 1]
    # resultsTheta[1, 2] <- 1/resultsAlpha[1, 1] * resultsAlpha[1, 2]
    colnames(resultsAlpha) <- c("Estimate", "Std. Error")
    rownames(resultsAlpha) <- "alpha"
    cat("\n   Overdispersion parameter alpha:\n\n")
    print(resultsAlpha)
    
    resultsVar <- matrix(resultsVar[-1, ], length(names(var.est0)[-1]), 4)
    rownames(resultsVar) <- names(var.est0)[-1]
    colnames(resultsVar)  <- c("Estimate", "Std. Error", "z value", "Pr(>z)")
    cat("\n\n   One sided Wald tests for variance components:\n\n")
    print(resultsVar)
  }
  
  tbr <- list(coefficients = list(fixed = coef0, random = ran.eff0), var.est = var.est0, std.err = c(std.err0, std.err1), z.val = c(zval0, zval1))
  invisible(tbr)
}