print.yuen <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nTest statistic: ", round(x$test, 4), " (df = ", round(x$df, 2), ")",", p-value = ", round(x$p.value, 5), "\n", sep = "")
  cat("\nTrimmed mean difference: ", round(x$diff, 5), "\n")
  cat("95 percent confidence interval:\n")
  cat(round(x$conf.int[1], 4), "   ", round(x$conf.int[2], 4), "\n\n")
}
