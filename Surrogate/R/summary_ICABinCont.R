summary.ICA.BinCont <- function(object, ..., Object){
 
  options(digits = 4)
  
  if (missing(Object)){Object <- object} 
  
  Object$R2_H <- na.exclude(Object$R2_H)

  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n# Total number of valid R2_H values")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(length(Object$R2_H))
  
  cat("\n\n\n# R2_H results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) R2_H: ", format(round(mean(Object$R2_H), 4), nsmall = 4), " (", format(round(sd(Object$R2_H), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$R2_H), 4), nsmall = 4), "; max: ",  format(round(max(Object$R2_H), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(Object$R2_H)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the R2_H distribution: \n\n")
  quant <- quantile(Object$R2_H, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
}
  

