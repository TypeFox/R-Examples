summary.FixedContContIT <- function(object, ..., Object){

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
  
  means_table <- rbind(tapply(Object$Data.Analyze$Surr, list(Object$Data.Analyze$Treat), mean), tapply(Object$Data.Analyze$True, list(Object$Data.Analyze$Treat), mean))
  colnames(means_table) <- c("Control Treatment", "Experimental treatment")
  rownames(means_table) <- c("Surrogate", "True endpoint")
  cat("\n\nMean surrogate and true endpoint values in each treatment group: \n\n")
  print(format(round(data.frame(means_table), 4), nsmall = 4))
  Var_table <- rbind(tapply(Object$Data.Analyze$Surr, list(Object$Data.Analyze$Treat), var), tapply(Object$Data.Analyze$True, list(Object$Data.Analyze$Treat), var))
  colnames(Var_table) <- c("Control Treatment", "Experimental treatment")
  rownames(Var_table) <- c("Surrogate", "True endpoint")
  cat("\n\nVar surrogate and true endpoint values in each treatment group: \n\n")
  print(format(round(data.frame(Var_table), 4), nsmall = 4))
  
  cat("\n\nCorrelations between the true and surrogate endpoints in the control (r_T0S0)")
  cat("\nand the experimental treatment groups (r_T1S1):\n\n")
  print(round(Object$Cor.Endpoints, 4), nsmall = 4)
  cat("\n\n\n# Information-theoretic surrogacy estimates summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  cat("Trial-level surrogacy (R2_ht): \n")
  print(format(round(Object$R2ht, 4), nsmall = 4))
  cat("\nIndividual-level surrogacy (R2_h.ind.clust): \n")
  print(format(round(Object$R2h.ind.clust, 4), nsmall = 4))
  cat("\nIndividual-level surrogacy assuming N=1 (R2_h.ind): \n")
  print(format(round(Object$R2h.ind, 4), nsmall = 4))
}



summary.MixedContContIT <- function(object, ..., Object){
  
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
  
  means_table <- rbind(tapply(Object$Data.Analyze$Surr, list(Object$Data.Analyze$Treat), mean), tapply(Object$Data.Analyze$True, list(Object$Data.Analyze$Treat), mean))
  colnames(means_table) <- c("Control Treatment", "Experimental treatment")
  rownames(means_table) <- c("Surrogate", "True endpoint")
  cat("\n\nMean surrogate and true endpoint values in each treatment group: \n\n")
  print(format(round(data.frame(means_table), 4), nsmall = 4))
  Var_table <- rbind(tapply(Object$Data.Analyze$Surr, list(Object$Data.Analyze$Treat), var), tapply(Object$Data.Analyze$True, list(Object$Data.Analyze$Treat), var))
  colnames(Var_table) <- c("Control Treatment", "Experimental treatment")
  rownames(Var_table) <- c("Surrogate", "True endpoint")
  cat("\n\nVar surrogate and true endpoint values in each treatment group: \n\n")
  print(format(round(data.frame(Var_table), 4), nsmall = 4))
  
  cat("\n\nCorrelations between the true and surrogate endpoints in the control (r_T0S0)")
  cat("\nand the experimental treatment groups (r_T1S1):\n\n")
  print(round(Object$Cor.Endpoints, 4), nsmall = 4)
  cat("\n\n\n# Information-theoretic surrogacy estimates summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  cat("Trial-level surrogacy (R2_ht): \n")
  print(format(round(Object$R2ht, 4), nsmall = 4))
  cat("\nIndividual-level surrogacy (R2_hind): \n")
  print(format(round(Object$R2h.ind, 4), nsmall = 4))
}


