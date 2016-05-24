summary.Explore.WS.Corr <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  x <- Object
  cat("Estimated correlations (as a function of time lag):") 
  cat("\n#==================================================\n")
  print(x$Est.Corr)
  cat("\n", (1-x$Alpha)*100, "% bootstrapped Confidence Intervals", sep="") 
  cat("\nLower limit:")
  cat("\n#=======================================\n")
  print(x$CI.Lower)
  cat("\nUpper limit:")
  cat("\n#=======================================\n")
  print(x$CI.Upper)
}
