summary.GoodPretreatContCont <- function(object, ..., Object){
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n\n# Rho2.Min results summary (Inf values are excluded)")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) Rho^2_min: ", format(round(mean(Object$Rho2.Min[which(is.finite(Object$Rho2.Min))]), 4), nsmall = 4), " (", format(round(sd(Object$Rho2.Min[which(is.finite(Object$Rho2.Min))]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$Rho2.Min[which(is.finite(Object$Rho2.Min))]), 4), nsmall = 4), "; max: ",  format(round(max(Object$Rho2.Min[which(is.finite(Object$Rho2.Min))]), 4), nsmall = 4), "]", sep="")
  cat("\n\nQuantiles of the Rho2.Min distribution: \n\n")
  quant <- quantile(Object$Rho2.Min, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)

  if (min(quant)<0){
    cat('\nNote. Some Rho2.Min values were negative. This indicates that the PMSE is so large that any Rho2.Min value')
    cat('\nsuffices to achieve the desired prediction accuracy.')
    
  }
  
}
