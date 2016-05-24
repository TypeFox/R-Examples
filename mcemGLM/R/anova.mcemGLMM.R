anova.mcemGLMM <- function(object, ...) {
#   if (!missing(opt)) {
#     if (class(opt) != "mcemGLMM") {
#       stop("Wrong class object.")
#     }
#     if (!(prod(colnames(opt$mcemEST) %in% colnames(object$mcemEST)) | prod(colnames(object$mcemEST) %in% colnames(opt$mcemEST)))) {
#       cat("    The models are not nested. Likelihood ratio test is not appropriate.\n\n")
#     }
#     
#     # Perform a likelihood ratio test.
#     value <- 2 * abs(tail(opt$loglikeVal, 1) - tail(object$loglikeVal, 1))
#     df <- abs(ncol(opt$mcemEST) - ncol(object$mcemEST))
#     pval <- 1 - pchisq(value, df)
#     
#     cat(paste("  Test statistic value:", value))
#     cat(paste("\n  Degrees of freedom:", df))
#     cat(paste("\n  p value:", round(pval, 8)))
#     
#     return(invisible(list(Chi.Sq = value, Df = df, p.value = pval)))
#   }
  # Fixed effects
  coef0 <- tail(object$mcemEST, n = 1)[1:ncol(object$x)]
  names(coef0) <- colnames(object$mcemEST)[1:ncol(object$x)]
    
  # Covariance matrix and standard errors
  cmat0 <- solve(object$iMatrix)
  # std.err0 <- sqrt(diag(cmat0))
  
  # z2 values
  
  # Coefficients positions
  pred0 <- attr(object$x, "assign")
  
  imat0 <- solve(object$iMatrix)
  if (attr(terms(as.formula(object$call$fixed)), "intercept") == 0) {
    names0 <- attr(terms(as.formula(object$call$fixed)), "term.labels")
  } else {
    names0 <- c("(Intercept)", attr(terms(as.formula(object$call$fixed)), "term.labels"))
  }
  
  # Drop intercept
  if(pred0[1] == 0) {
    coef0 <- coef0[-1]
    pred0 <- pred0[-1]
    imat0 <- imat0[-1, -1]
    names0 <- names0[-1]
  }
  
  # Number of predictors 
  npred <- unique(pred0)
  wald0 <- rep(0, length(npred))
  for (i in npred) {
    ind0 <- which(pred0 == i)
    wald0[i] <- t(coef0[ind0]) %*% solve(imat0[ind0, ind0]) %*% coef0[ind0]
  }
  
  df0 <- table(pred0)
  pval0 <- round(pchisq(wald0, df0, lower.tail = FALSE), 8)
  tbr <- cbind(df0, wald0, pval0)
  cat("   Wald's Chi-squared ANOVA table\n\n")
  colnames(tbr) <- c("Df", "Wald Stat", "Pr(>W)")
  rownames(tbr) <- names0
  return(tbr)
}