print.med1way <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nTest statistic:", round(x$test, 4),"\n")
  cat("Critical value:", round(x$crit.val, 4),"\n")
  cat("p-value:", round(x$p.value, 5), "\n")
  cat("\n")
}
