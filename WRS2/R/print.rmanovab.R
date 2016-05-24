print.rmanovab <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nTest statistic:", round(x$test, 4),"\n")
  cat("Critical value:", round(x$crit, 4), "\n")
  if (x$test < x$crit) sig <- FALSE else sig <- TRUE
  cat("Significant: ", sig, "\n")
  cat("\n")
}
