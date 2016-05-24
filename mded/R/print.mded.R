print.mded <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
  distr.cases <- c(length(x$distr1), length(x$distr2))
  means.df <- data.frame(means = x$means, n = as.character(distr.cases))
  rownames(means.df) <- x$distr.names

  n <- matrix(c(x$cases, sum(x$cases)), nrow = 3)
  dimnames(n) <- list(c("true", "false", "total"), c("n"))

  cat("\nTest:\n")
  cat("H0 ", x$distr.names[1], "=", x$distr.names[2], "\n")
  cat("H1 ", x$distr.names[1], ">", x$distr.names[2], "\n")
  cat("significance level =", format.pval(x$stat, digits = digits), "\n")

  cat("\nData:\n")
  cat("distr1 =", x$distr.names[1], "\n")
  cat("distr2 =", x$distr.names[2], "\n")

  cat("\nMeans:\n")
  print(means.df, digits = digits, print.gap = 2, ...)

  cat("\nCases in the difference:\n")
  print(n, print.gap = 2, ...)
  cat("\n")

  invisible(x)
}
