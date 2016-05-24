print.spp <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nTest statistic:", round(x$test, 4),"\n")
  cat("p-value:", round(x$p.value, 5), "\n")
  cat("\n")
}
