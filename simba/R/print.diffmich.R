"print.diffmich" <-
  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  nperm <- x$permutations  
  cat("\n")
  cat("Is difference in the Parameters significant?", "\n")
  cat("Significance is based on", nperm, "permutations", "\n")
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("Difference in a: ")
  cat(formatC(x$diffa, digits = digits), "\n")
  if (nperm) {
    cat("Significance:", format.pval(x$signifa, eps = 1/nperm), 
        "\n\n")
    outa <- quantile(x$permsa, c(0.9, 0.95, 0.975, 0.99))
    cat("Empirical upper confidence limits of a:\n")
    print(outa, digits = 3)
  }
  cat("\n\n")
  cat("Difference in b: ")
  cat(formatC(x$diffb, digits = digits), "\n")
  if (nperm) {
    cat("Significance:", format.pval(x$signifb, eps = 1/nperm), 
        "\n\n")
    outb <- quantile(x$permsb, c(0.9, 0.95, 0.975, 0.99))
    cat("Empirical upper confidence limits of r:\n")
    print(outb, digits = 3)
  }
  cat("\n\n")
  invisible(x)
}

