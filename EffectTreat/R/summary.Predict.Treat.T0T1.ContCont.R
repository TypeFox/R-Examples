summary.Predict.Treat.T0T1.ContCont <- function(object, ..., Object){
  if (missing(Object)){Object <- object} 
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Predicted (Mean) Delta_T_j | S_j")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(Object$Pred_T)
  CI_Val <- 1-Object$alpha
  cat("\n\n\n# Variance and ", CI_Val*100, "% support interval of Delta_T_j | S_j for the requested value of rho_T0T1", sep="")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  rho_T0T1 <- Object$T0T1
  Var_Delta.T <- Object$Var_Delta.T
  cat("                  rho_T0T1  ", "Var Delta_T_j | S_j  ", CI_Val*100, "% supp. int. around Delta_T_j | S_j \n\n")
  
  cat("Requested value:     ", format(round(rho_T0T1, 3), nsmall = 3), "            ", 
      format(round(Object$Var_Delta.T_S, 3), nsmall = 3), "          [",
      Object$CI_low, "; ", Object$CI_high, "]\n", sep="")
  
  aantal_Zero_included <- length(Object$Var_Delta.T[Object$CI_low<0 & Object$CI_high>0])
  aantal_CI_under_Zero <- length(Object$Var_Delta.T[Object$CI_low<0 & Object$CI_high<0])
  aantal_CI_boven_Zero <- length(Object$Var_Delta.T[Object$CI_low>0 & Object$CI_high>0])
  
  if (aantal_Zero_included>0){
    cat("\n => 0 is included in the support interval.\n")
  }
  
  if (aantal_CI_under_Zero>0){
    cat("\n => Entire support interval falls below 0. ")
  }
  if (aantal_CI_boven_Zero > 0){
    cat("\n => Entire support interval falls above 0. ")
  }
  
}


