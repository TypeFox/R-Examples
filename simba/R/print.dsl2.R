"print.dsl2" <-
  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  nperm <- x$permutations  
  cat("\n")
  cat("Is difference in slope significant?", "\n")
  cat("Significance is based on", nperm, "permutations", "\n")
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("Difference in Slope: ")
  cat(formatC(x$slope.diff, digits = digits), "\n")
  if (nperm) {
    cat("Significance:", format.pval(x$signif, eps = 1/nperm), 
        "\n")
    out <- quantile(x$perms, c(0.9, 0.95, 0.975, 0.99))
    cat("Empirical upper confidence limits of r:\n")
    print(out, digits = 3)
  }
  cat("\n")
  cat("Difference in Intercept: ")
  cat(formatC(x$intercept, digits = digits), "\n")
  if (nperm) {
    cat("Significance:", format.pval(x$signific, eps = 1/nperm), 
        "\n")
    out <- quantile(x$permsic, c(0.9, 0.95, 0.975, 0.99))
    cat("Empirical upper confidence limits of r:\n")
    print(out, digits = 3)
  }
  cat("\n")
  invisible(x)
}

