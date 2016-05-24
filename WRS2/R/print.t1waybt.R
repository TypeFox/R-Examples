print.t1waybt <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nEffective number of bootstrap samples was ", x$nboot.eff, ".\n", sep = "")
  cat("\nTest statistic:", round(x$test, 4),"\n")
  cat("p-value:", round(x$p.value, 5), "\n")
  cat("Variance explained", round(x$Var.Explained, 3), "\n")
  cat("Effect size", round(x$Effect.Size, 3), "\n")
  cat("\n")
}
