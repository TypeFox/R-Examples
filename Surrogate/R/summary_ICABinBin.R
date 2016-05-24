summary.ICA.BinBin <- function(object, ..., Object){
 
  options(digits = 4)
  
  if (missing(Object)){Object <- object} 
  
  if (Object$Monotonicity!="General"){
  
#  if (missing(Object)){Object <- object} 
  
  Object$R2_H <- na.exclude(Object$R2_H)
  Object$Theta_T <- na.exclude(Object$Theta_T)
  Object$Theta_S <- na.exclude(Object$Theta_S)
  
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n# Total number of valid Pi vectors")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(dim(Object$Pi.Vectors)[1])
  
  cat("\n\n\n# R2_H results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) R2_H: ", format(round(mean(Object$R2_H), 4), nsmall = 4), " (", format(round(sd(Object$R2_H), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$R2_H), 4), nsmall = 4), "; max: ",  format(round(max(Object$R2_H), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(Object$R2_H)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the R2_H distribution: \n\n")
  quant <- quantile(Object$R2_H, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  cat("\n\n\n# R_H results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) R_H: ", format(round(mean(sqrt(Object$R2_H)), 4), nsmall = 4), " (", format(round(sd(sqrt(Object$R2_H)), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(sqrt(Object$R2_H)), 4), nsmall = 4), "; max: ",  format(round(max(sqrt(Object$R2_H)), 4), nsmall = 4), "]", sep="")
  cat("\nMode R_H: ", format(round(mode(sqrt(Object$R2_H))$mode_val, 4), nsmall = 4))
  
  cat("\n\nQuantiles of the R_H distribution: \n\n")
  quant <- quantile(sqrt(Object$R2_H), probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  
#  cat("\n\n# C3 results summary")
#  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
#  cat("Mean (SD) C3: ", format(round(mean(Object$C3), 4), nsmall = 4), " (", format(round(sd(Object$C3), 4), nsmall = 4), ")", 
#      "  [min: ", format(round(min(Object$C3), 4), nsmall = 4), "; max: ",  format(round(max(Object$C3), 4), nsmall = 4), "]", sep="")
#  cat("\nMode C3: ", format(round(mode(Object$C3)$mode_val, 4), nsmall = 4))
#  cat("\n\nQuantiles of the C3 distribution: \n\n")
#  quant <- quantile(Object$C3, probs = c(.05, .10, .20, .50, .80, .90, .95))
#  print(quant)
  
  cat("\n\n# Theta_T results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) Theta_T: ", format(round(mean(Object$Theta_T), 4), nsmall = 4), " (", format(round(sd(Object$Theta_T), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$Theta_T), 4), nsmall = 4), "; max: ",  format(round(max(Object$Theta_T), 4), nsmall = 4), "]", sep="")
  if (Object$Theta_T[1]!=Inf) {
    cat("\nMode Theta_T: ", format(round(mode(Object$Theta_T)$mode_val, 4), nsmall = 4))}
  cat("\n\nQuantiles of the Theta_T distribution: \n\n")
  quant <- quantile(Object$Theta_T, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  cat("\n\n# Theta_S results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) Theta_S: ", format(round(mean(Object$Theta_S), 4), nsmall = 4), " (", format(round(sd(Object$Theta_S), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$Theta_S), 4), nsmall = 4), "; max: ",  format(round(max(Object$Theta_S), 4), nsmall = 4), "]", sep="")
  if (Object$Theta_S[1]!=Inf) {
    cat("\nMode Theta_S: ", format(round(mode(Object$Theta_S)$mode_val, 4), nsmall = 4))}
  cat("\n\nQuantiles of the Theta_S distribution: \n\n")
  quant <- quantile(Object$Theta_S, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  }

if (Object$Monotonicity=="General"){

  Object$R2_H <- na.exclude(Object$R2_H)
  Object$Theta_T <- na.exclude(Object$Theta_T)
  Object$Theta_S <- na.exclude(Object$Theta_S)
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n# Number of valid Pi vectors")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat("Total:", dim(Object$Pi.Vectors)[1])
  cat("\n\nIn the different montonicity scenarios:")
  print(table(Object$Pi.Vectors$Monotonicity))
  
  results <- cbind.data.frame(Object$Pi.Vectors$Monotonicity, Object$R2_H)
  
  cat("\n\n# Summary of results obtained in different monotonicity scenarios")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  
  cat("# R2_H results summary")
  cat("\n~~~~~~~~~~~~~~~~~~~~~\n\nMean:\n")
  print(tapply(results[,2], list(results[,1]), mean))
  
  pc50 <- function(x=x){
    quantile(x = x, probs = .5, na.rm = T)}
  cat("\nMedian:\n")
  print(tapply(results[,2], list(results[,1]), pc50))
  
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  cat("\nMode:\n")
  print(data.frame(t(tapply(results[,2], list(results[,1]), mode)), row.names=""))
  
  cat("\nSD:\n")
  print(tapply(results[,2], list(results[,1]), sd))
  
  cat("\nMin:\n")
  print(tapply(results[,2], list(results[,1]), min))
  cat("\nMax:\n")
  print(tapply(results[,2], list(results[,1]), max))
  
  #  results <- cbind.data.frame(Object$Pi.Vectors$Monotonicity, Object$C3)
  
  #  cat("\n\n# C3 results summary")
  #  cat("\n#~~~~~~~~~~~~~~~~~~~\n\nMean:\n")
  #  print(tapply(results[,2], list(results[,1]), mean))
  
  #  cat("\nMedian:\n")
  #  print(tapply(results[,2], list(results[,1]), pc50))
  
  #  cat("\nMode:\n")
  #  print(data.frame(t(tapply(results[,2], list(results[,1]), mode)), row.names=""))
  
  #  cat("\nSD:\n")
  #  print(tapply(results[,2], list(results[,1]), sd))
  
  #  cat("\nMin:\n")
  #  print(tapply(results[,2], list(results[,1]), min))
  #  cat("\nMax:\n")
  #  print(tapply(results[,2], list(results[,1]), max))
  
  
  results <- cbind.data.frame(Object$Pi.Vectors$Monotonicity, Object$Theta_T)
  
  cat("\n\n# Theta_T results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~\n\nMean:\n")
  print(tapply(results[,2], list(results[,1]), mean))
  
  cat("\nMedian:\n")
  print(tapply(results[,2], list(results[,1]), pc50))
  
  cat("\nSD:\n")
  print(tapply(results[,2], list(results[,1]), sd))
  
  cat("\nMin:\n")
  print(tapply(results[,2], list(results[,1]), min))
  cat("\nMax:\n")
  print(tapply(results[,2], list(results[,1]), max))
  
  
  results <- cbind.data.frame(Object$Pi.Vectors$Monotonicity, Object$Theta_S)
  
  cat("\n\n# Theta_S results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~\n\nMean:\n")
  print(tapply(results[,2], list(results[,1]), mean))
  
  cat("\nMedian:\n")
  print(tapply(results[,2], list(results[,1]), pc50))
  
  cat("\nSD:\n")
  print(tapply(results[,2], list(results[,1]), sd))
  
  cat("\nMin:\n")
  print(tapply(results[,2], list(results[,1]), min))
  cat("\nMax:\n")
  print(tapply(results[,2], list(results[,1]), max))
  

}
}
  

