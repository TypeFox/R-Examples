"print.pclist" <-
function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  nperm <- x$permutations  
  cat("\n")
  cat("Correlation based on", x$method, "\n")
  cat("Mantel-like approach, based on", nperm, "permutations", "\n")
  cat("on each level", "\n")
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("Overall number of cases:", x$gesN, "\n\n")
  cat("Table of permuted correlation on", x$strata, "levels:", "\n")
  print(x$out, digits = 4)
  cat("\n\n")
  invisible(x)
}