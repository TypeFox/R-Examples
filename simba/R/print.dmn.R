"print.dmn" <-
  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  nperm <- x$permutations  
  cat("\n")
  cat("Are difference in Mean and F significant?", "\n")
  cat("Significance is based on", nperm, "permutations", "\n")
  cat("\nCall: ")
  cat(deparse(x$call), "\n\n")
  cat("Overall Mean: ")
  cat(formatC(x$meanom, digits = digits), "\n\n")
  cat("mean(x): ")
  cat(formatC(x$mean.x, digits = digits), "\n")
  cat("mean(y): ")
  cat(formatC(x$mean.y, digits = digits), "\n")
  cat("Difference in Mean: ")
  cat(formatC(x$diff, digits = digits), "\n\n")
  if (nperm) {
    cat("Significance of Difference of Means:", format.pval(x$sig, eps = 1/nperm), 
        "\n\n")
    outM <- quantile(x$bootsM, c(0.9, 0.95, 0.975, 0.99))
    cat("Empirical upper confidence limits of Means:\n")
    print(outM, digits = 3)
  }
  cat("\n\n")
  cat("F-value: ")
  cat(formatC(x$F, digits=3), "\n\n")
  if (nperm) {
    cat("Significance of F:", format.pval(x$sigF, eps = 1/nperm), 
        "\n\n")
    outF <- quantile(x$bootsF, c(0.9, 0.95, 0.975, 0.99))
    cat("Empirical upper confidence limits of Fs:\n")
    print(outF, digits = 3)
  }
  cat("\n\n")
  invisible(x)
}

