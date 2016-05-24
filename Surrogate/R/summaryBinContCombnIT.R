summary.FixedBinBinIT <- summary.FixedBinContIT <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", nrow(Object$Obs.Per.Trial))
  cat("\nTotal number of patients: ", dim(Object$Data.Analyze)[1])
  cat("\nM(SD) patients per trial: ", format(round(mean((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), " (", format(round(sd((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), ")", 
      "  [min: ", min((Object$Obs.Per.Trial$Obs.per.trial)), "; max: ",  max((Object$Obs.Per.Trial$Obs.per.trial)), "]", sep="")
  cat("\nTotal number of patients in experimental treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat==1]), 
      "\nTotal number of patients in control treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat!=1])) 
  
  cat("\n\n\n# Information-theoretic surrogacy estimates summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  cat("Trial-level surrogacy (R2_ht): \n")
  print(format(round(Object$R2ht, 4), nsmall = 4))
#  cat("\nIndividual-level surrogacy (R2_h): \n")
#  print(format(round(Object$R2h, 4), nsmall = 4))
  cat("\nIndividual-level surrogacy (R2_h.ind): \n")
  print(format(round(Object$R2h.ind, 4), nsmall = 4))
#  cat("\nIndividual-level surrogacy taking max. bound into account (R2_b.ind): \n")
#  print(format(round(Object$R2b.ind, 4), nsmall = 4))
}



summary.FixedContBinIT <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", nrow(Object$Obs.Per.Trial))
  cat("\nTotal number of patients: ", dim(Object$Data.Analyze)[1])
  cat("\nM(SD) patients per trial: ", format(round(mean((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), " (", format(round(sd((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), ")", 
      "  [min: ", min((Object$Obs.Per.Trial$Obs.per.trial)), "; max: ",  max((Object$Obs.Per.Trial$Obs.per.trial)), "]", sep="")
  cat("\nTotal number of patients in experimental treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat==1]), 
      "\nTotal number of patients in control treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat!=1])) 
  
  cat("\n\n\n# Information-theoretic surrogacy estimates summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  cat("Trial-level surrogacy (R2_ht): \n")
  print(format(round(Object$R2ht, 4), nsmall = 4))
  cat("\nIndividual-level surrogacy (R2_h.ind): \n")
  print(format(round(Object$R2h.ind, 4), nsmall = 4))
}