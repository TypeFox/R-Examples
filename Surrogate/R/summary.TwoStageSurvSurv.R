
summary.TwoStageSurvSurv <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", length(Object$Results.Stage.1$Trial.Name))
  cat("\nM(SD) trial size: ", format(round(mean((Object$Results.Stage.1$Trial.Size)), 4), nsmall = 4),
      " (", format(round(sd((Object$Results.Stage.1$Trial.Size)), 4), nsmall = 4), ")", 
      "  [min: ", min((Object$Results.Stage.1$Trial.Size)), "; max: ",  max((Object$Results.Stage.1$Trial.Size)), "]", sep="")
  
  cat("\n\n\n# R^2_{ht} results")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  print(format(round(Object$R2.ht, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$R.ht, 4), nsmall = 4))
  cat("\n")
  }
