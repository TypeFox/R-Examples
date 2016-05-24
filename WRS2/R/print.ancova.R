print.ancova <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  if (names(x)[5] == "se") {
   df <- data.frame(x$n1, x$n2, x$trDiff, x$se, x$ci.low, x$ci.hi, x$test,  x$p.vals)
   rownames(df) <- paste(x$cnames[3], "=", x$evalpts)
   colnames(df) <- c(paste("n:", x$faclevels[1]), paste("n:",x$faclevels[2]), "trimmed mean diff", "se", "lower CI", "upper CI", "statistic", "p-value")
  } else {
    df <- data.frame(x$n1, x$n2, x$trDiff, x$ci.low, x$ci.hi, x$test,  x$p.vals)
    rownames(df) <- paste(x$cnames[3], "=", x$evalpts)
    colnames(df) <- c(paste("n:", x$faclevels[1]), paste("n:",x$faclevels[2]), "trimmed mean diff", "lower CI", "upper CI", "statistic", "p-value")
  } 
  print(round(df, 4)) 
  cat("\n")
}
