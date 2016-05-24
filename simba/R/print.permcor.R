"print.permcor" <-
function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  nperm <- x$permutations  
  cat("\n")
  cat("Correlation based on", x$method, "\n")
  cat("Mantel-like approach, based on", nperm, "permutations", "\n")
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("Correlation r: ")
  cat(formatC(x$statistic, digits = digits), " ( n=", x$n, ")", "\n")
  if (nperm) {
    cat("Significance:", format.pval(x$signif, eps = 1/nperm), 
        "\n\n")
    out <- quantile(x$perms, c(0.9, 0.95, 0.975, 0.99))
    cat("Empirical upper confidence limits of r:\n")
    print(out, digits = 3)
  }
  cat("\n\n")
  invisible(x)
}

