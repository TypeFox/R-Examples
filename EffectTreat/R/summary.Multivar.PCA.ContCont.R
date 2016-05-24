summary.Multivar.PCA.ContCont <- function(object, ..., Object){
  if (missing(Object)){Object <- object} 
  
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Total number of matrices that can be formed by the specified vectors and/or scalars")
  cat("\n# of the correlations in the function call")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(Object$Total.Num.Matrices)
  cat("\n\n# Total number of positive definite matrices")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(Object$Pos.Def)
  cat("\n\n\n# Predictive causal association (PCA; R^2_{psi}) results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) PCA: ", format(round(mean(Object$PCA), 4), nsmall = 4), " (", format(round(sd(Object$PCA), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$PCA), 4), nsmall = 4), "; max: ",  format(round(max(Object$PCA), 4), nsmall = 4), "]", sep="")
  cat("\nMode PCA: ", format(round(mode(Object$PCA)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the PCA distribution: \n\n")
  quant <- quantile(Object$PCA, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
}
