summary.SurvSurv <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
#  cat("\n\n# Data summary and descriptives")
#  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#  cat("\nM(SD) trial size: ", format(round(mean((Object$Results.Stage.1$Trial.Size)), 4), nsmall = 4),
#      " (", format(round(sd((Object$Results.Stage.1$Trial.Size)), 4), nsmall = 4), ")", 
#      "  [min: ", min((Object$Results.Stage.1$Trial.Size)), "; max: ",  max((Object$Results.Stage.1$Trial.Size)), "]", sep="")
  
  cat("\n\n\n# R^2_{trial} results")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n")
  print(format(round(Object$R2.trial, 4), nsmall = 4))
  cat("\n")
  
  cat("\n\n# R^2_{h.ind} (LRF) results")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\nOverall R^2_{h.ind} (LRF): \n\n")
  print(format(round(Object$R2.hind, 4), nsmall = 4))
  
  cat("\n\n# R^2_{h.ind.QF} (LRF_a; O'Quinly and Flandre, 2006) results")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\nOverall R^2_{h.ind.QF} (LRF_a): \n\n")
  print(format(round(Object$R2h.ind.QF, 4), nsmall = 4))
  
  
  cat("\n\nR^2_{h.ind.QF} (LRF_a) per trial: \n\n")
  print(format(round(Object$R2.hInd.By.Trial.QF, 4), nsmall = 4))
  cat("\n")
  }
