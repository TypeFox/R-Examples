print.pb2 <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nTest statistic: ", round(x$test, 4), ", p-value = ", round(x$p.value, 5), "\n", sep = "")
  cat("95 percent confidence interval:\n")
  cat(round(x$conf.int[1], 4), "   ", round(x$conf.int[2], 4), "\n\n")
}
