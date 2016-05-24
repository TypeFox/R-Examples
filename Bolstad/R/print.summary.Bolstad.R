#' @export

print.summary.Bolstad = function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n")
  cat(paste0(deparse(x$call), sep = "\n"), "\n")
  
  cat("Residuals:", "\n")
  resid = x$residuals
  if (x$res.df > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(resid)) == 2L) 
      structure(apply(t(resid), 1L, quantile), dimnames = list(nam, 
                                                               dimnames(resid)[[2L]]))
    else {
      zz <- zapsmall(quantile(resid), digits + 1L)
      structure(zz, names = nam)
    }
  }
  print(rq, digits = digits, ...)
  
  cat("\nCoefficients:\n")
  coef.mat = matrix(nrow = x$rank, ncol = 3)
  
  colnames(coef.mat) = c("Posterior Mean", "Std. Error", "t value")
  rownames(coef.mat) = c(x$terms)
  
  coef.mat[, "Posterior Mean"] = x$coef
  coef.mat[, "Std. Error"] = x$std.err
  coef.mat[, "t value"] = x$coef / x$std.err
  
  print(format(coef.mat, digits = digits), print.gap = 2L, quote = FALSE, ...)
  
  cat("---\n")
  
  if (!is.null(x$prior)) {
    prior.coef = x$prior$b0
    cat("Prior Coefficients:\n")
    names(prior.coef) = x$terms
    print(format(prior.coef, digits = digits), print.gap = 2L, quote = FALSE, ...)

    prior.cov = x$prior$V0
    cat("\nPrior Covariance Matrix:\n")
    dimnames(prior.cov) = list(x$terms, x$terms)
    print(format(prior.cov, digits = digits), print.gap = 2L, quote = FALSE, ...)
  } 
  else {
    cat("\nNote: No prior given (Using flat prior).\n")
  }
}