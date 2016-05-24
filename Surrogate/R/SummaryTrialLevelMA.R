summary.TrialLevelMA <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", length(Object$N.Vector))

  cat("\n\n\n# Meta-analytic results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  print(format(round(Object$Trial.R2, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Trial.R, 4), nsmall = 4))
}

