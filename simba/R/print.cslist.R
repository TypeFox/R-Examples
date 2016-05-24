"print.cslist" <-
  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  nperm <- x$permutations  
  cat("\n")
  cat("Differences in mean based on", x$method, "\n")
  cat("Inference obtained with", nperm, "permutations", "\n")
  cat("for each level-comparison", "\n")
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("Overall number of cases:", x$gesN, "\n\n")
  cat("Table of permuted correlation on", x$strata, "levels:", "\n")
  print(x$out, digits = 4)
  cat("\n\n")
  cat("Table of Connections", "\n")
  print(x$cons)
  cat("\n")
  invisible(x)
}