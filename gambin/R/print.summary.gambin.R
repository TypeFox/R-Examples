print.summary.gambin <-
function(x, ...) 
{
  cat("\n")
  cat("GamBin distribution fit\n")
  cat("Data: ")
  cat(x$Dataname)
  cat("\n\n")
  
  vals <- data.frame(Estimated = x$Alpha, CI95_low = x$ConfInt95[1], CI95_high = x$ConfInt95[2])
  rownames(vals) <- "Alpha"
  
  cat("Coefficients:\n")
  print(vals)

  cat("MaxOctave\t")
  cat(floor(x$MaxOctave))
  
  cat("\n\n")
  cat(paste("Chi-squared fit:", paste("X^2 =", round(x$ChiSq$statistic,3)), paste ("df =", x$Fit$parameter), paste("p-value =", round(x$ChiSq$p.value,3)), "\n", sep = "\t" ))
  cat("\n")
}
