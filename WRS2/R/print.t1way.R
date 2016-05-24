print.t1way <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nTest statistic:", round(x$test, 4),"\n")
  cat("Degrees of Freedom 1:", round(x$df1, 2),"\n")
  cat("Degrees of Freedom 2:", round(x$df2, 2) ,"\n")
  cat("p-value:", round(x$p.value, 5), "\n")
  cat("\n")
}
