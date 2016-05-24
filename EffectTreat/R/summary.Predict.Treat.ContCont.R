summary.Predict.Treat.ContCont <- summary.Predict.Treat.Multivar.ContCont <- function(object, ..., Object){
  if (missing(Object)){Object <- object} 
  
  if (class(Object)=="Predict.Treat.ContCont"){
    
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Predicted (Mean) Delta_T_j | S_j")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(Object$Pred_T)
  cat("\n\n\n# Variances and 95% support intervals of Delta_T_j | S_j for different values of rho_T0T1")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  min_rho_T0T1 <- min(Object$T0T1)
  max_rho_T0T1 <- max(Object$T0T1)
  median_rho_T0T1 <- median(Object$T0T1)
  mean_rho_T0T1 <- mean(Object$T0T1)
  
  min_Var_Delta.T <- Object$Var_Delta.T[Object$T0T1==min_rho_T0T1]
  max_Var_Delta.T <- Object$Var_Delta.T[Object$T0T1==max_rho_T0T1]
  median_Var_Delta.T <- Object$Var_Delta.T[which.min(abs(Object$T0T1-median_rho_T0T1))]
  mean_Var_Delta.T <- Object$Var_Delta.T[which.min(abs(Object$T0T1-mean_rho_T0T1))]
  
  min_Var_Delta.T_givenS <- min_Var_Delta.T*(1-((Object$PCA[Object$T0T1==min_rho_T0T1])**2))
  max_Var_Delta.T_givenS <- max_Var_Delta.T*(1-((Object$PCA[Object$T0T1==max_rho_T0T1])**2))
  median_Var_Delta.T_givenS <- median_Var_Delta.T*(1-((Object$PCA[which.min(abs(Object$T0T1-median_rho_T0T1))])**2))
  mean_Var_Delta.T_givenS <- mean_Var_Delta.T*(1-((Object$PCA[which.min(abs(Object$T0T1-mean_rho_T0T1))])**2))
  rm(min_Var_Delta.T, max_Var_Delta.T, median_Var_Delta.T, mean_Var_Delta.T)
  
  crit_val <- qnorm(c(.05/2), mean=0, sd=1, lower.tail=FALSE)
  
  cat("                  rho_T0T1  ", "Var Delta_T_j | S_j  ", "95% supp. int. around Delta_T_j | S_j \n\n")
  
  cat(" (min. value)     ", format(round(min_rho_T0T1, 3), nsmall = 3), "            ", 
      format(round(min_Var_Delta.T_givenS, 3), nsmall = 3), "          [",
      (-crit_val*sqrt(min_Var_Delta.T_givenS))+Object$Pred_T, "; ", (crit_val*sqrt(min_Var_Delta.T_givenS))+Object$Pred_T, "]\n",
      sep="")
  
  cat(" (max. value)     ", format(round(max_rho_T0T1, 3), nsmall = 3), "            ", 
      format(round(max_Var_Delta.T_givenS, 3), nsmall = 3), "          [",
      (-crit_val*sqrt(max_Var_Delta.T_givenS))+Object$Pred_T, "; ", (crit_val*sqrt(max_Var_Delta.T_givenS))+Object$Pred_T, "]\n",
      sep="")
  
  cat(" (median value)   ", format(round(median_rho_T0T1, 3), nsmall = 3), "            ", 
      format(round(median_Var_Delta.T_givenS, 3), nsmall = 3), "          [",
      (-crit_val*sqrt(median_Var_Delta.T_givenS))+Object$Pred_T, "; ", (crit_val*sqrt(median_Var_Delta.T_givenS))+Object$Pred_T, "]\n",
      sep="")
  
  cat(" (mean value)     ", format(round(mean_rho_T0T1, 3), nsmall = 3), "            ", 
      format(round(mean_Var_Delta.T_givenS, 3), nsmall = 3), "          [",
      (-crit_val*sqrt(mean_Var_Delta.T_givenS))+Object$Pred_T, "; ", (crit_val*sqrt(mean_Var_Delta.T_givenS))+Object$Pred_T, "]\n",
      sep="")
  
  cat("\n\n\n# Proportion of 95% support intervals for Delta_T_j | S_j\n") 
  cat("that include 0, are < 0, and are > 0")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  
  CIs_LB <- (-crit_val*sqrt(Object$Var_Delta.T * (1-(Object$PCA**2))))+Object$Pred_T 
  CIs_UB <- (crit_val*sqrt(Object$Var_Delta.T * (1-(Object$PCA**2))))+Object$Pred_T 
  
  aantal_Zero_included <- length(Object$Var_Delta.T[CIs_LB<0 & CIs_UB>0])
  aantal_CI_under_Zero <- length(Object$Var_Delta.T[CIs_LB<0 & CIs_UB<0])
  aantal_CI_boven_Zero <- length(Object$Var_Delta.T[CIs_LB>0 & CIs_UB>0])
  
  temp_vect <- rep(NA, length(Object$Var_Delta.T))
  temp_vect[CIs_LB<0 & CIs_UB>0] <- 1 #zero included
  temp_vect[CIs_LB<0 & CIs_UB<0] <- 2  # under zero
  temp_vect[CIs_LB>0 & CIs_UB>0] <- 3  # above zero
  
  if (aantal_Zero_included>0){
    cat("0 included in support interval: ", aantal_Zero_included/length(Object$Var_Delta.T),
        "     (obtained for rho_T0T1 values in range [", min(Object$T0T1[temp_vect==1]), "; ", 
        max(Object$T0T1[temp_vect==1]), "])\n", sep="")
  }
  if (aantal_Zero_included==0){
    cat("0 included in support interval: ", aantal_Zero_included/length(Object$Var_Delta.T), "\n")
  }  
  
  if (aantal_CI_under_Zero>0){
    cat("Entire support interval below 0: ", aantal_CI_under_Zero/length(Object$Var_Delta.T),
        "     (obtained for rho_T0T1 values in range [", min(Object$T0T1[temp_vect==2]), "; ", 
        max(Object$T0T1[temp_vect==2]), "])\n", sep="")
  }
  if (aantal_CI_under_Zero==0){
    cat("Entire support interval below 0: ", aantal_CI_under_Zero/length(Object$Var_Delta.T), "\n")
  }
  
  if (aantal_CI_boven_Zero > 0){
    cat("Entire support interval above 0: ", aantal_CI_boven_Zero/length(Object$Var_Delta.T),
        "     (obtained for rho_T0T1 values in range [", min(Object$T0T1[temp_vect==3]), "; ", 
        max(Object$T0T1[temp_vect==3]), "])\n", sep="")
  }
  if (aantal_CI_boven_Zero == 0){
    cat("Entire support interval above 0: ", aantal_CI_boven_Zero/length(Object$Var_Delta.T), "\n")
  }
  }



  if (class(Object)=="Predict.Treat.Multivar.ContCont"){
    
    cat("\nFunction call:\n\n")
    print(Object$Call)
    cat("\n\n# Predicted (Mean) Delta_T_j | S_j")
    cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    cat(Object$Pred_T)
    cat("\n\n\n# Variances and 95% support intervals of Delta_T_j | S_j for different values of rho_T0T1")
    cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
    min_rho_T0T1 <- min(Object$T0T1)
    max_rho_T0T1 <- max(Object$T0T1)
    median_rho_T0T1 <- median(Object$T0T1)
    mean_rho_T0T1 <- mean(Object$T0T1)
    
    tab <- cbind(Object$T0T1, Object$Var_Delta.T_S)
    min_Var_Delta.T_givenS <- max(tab[,2][order(tab[,1])]) #LV
    max_Var_Delta.T_givenS <- min(tab[,2][order(tab[,1])]) #LV
    median_Var_Delta.T_givenS <- median(Object$Var_Delta.T_S)
    mean_Var_Delta.T_givenS <- mean(Object$Var_Delta.T_S)
  
    crit_val <- qnorm(c(.05/2), mean=0, sd=1, lower.tail=FALSE)
    
    cat("                  rho_T0T1  ", "Var Delta_T_j | S_j  ", "95% supp. int. around Delta_T_j | S_j \n\n")
    
    cat(" (min. value)     ", format(round(min_rho_T0T1, 3), nsmall = 3), "            ", 
        format(round(min_Var_Delta.T_givenS, 3), nsmall = 3), "          [",
        (-crit_val*sqrt(min_Var_Delta.T_givenS))+Object$Pred_T, "; ", (crit_val*sqrt(min_Var_Delta.T_givenS))+Object$Pred_T, "]\n",
        sep="")
    
    cat(" (max. value)     ", format(round(max_rho_T0T1, 3), nsmall = 3), "            ", 
        format(round(max_Var_Delta.T_givenS, 3), nsmall = 3), "          [",
        (-crit_val*sqrt(max_Var_Delta.T_givenS))+Object$Pred_T, "; ", (crit_val*sqrt(max_Var_Delta.T_givenS))+Object$Pred_T, "]\n",
        sep="")
    
    cat(" (median value)   ", format(round(median_rho_T0T1, 3), nsmall = 3), "            ", 
        format(round(median_Var_Delta.T_givenS, 3), nsmall = 3), "          [",
        (-crit_val*sqrt(median_Var_Delta.T_givenS))+Object$Pred_T, "; ", (crit_val*sqrt(median_Var_Delta.T_givenS))+Object$Pred_T, "]\n",
        sep="")
    
    cat(" (mean value)     ", format(round(mean_rho_T0T1, 3), nsmall = 3), "            ", 
        format(round(mean_Var_Delta.T_givenS, 3), nsmall = 3), "          [",
        (-crit_val*sqrt(mean_Var_Delta.T_givenS))+Object$Pred_T, "; ", (crit_val*sqrt(mean_Var_Delta.T_givenS))+Object$Pred_T, "]\n",
        sep="")
    
    cat("\n\n\n# Proportion of 95% support intervals for Delta_T_j | S_j\n") 
    cat("that include 0, are < 0, and are > 0")
    cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
    
    CIs_LB <- (-crit_val*sqrt(Object$Var_Delta.T_S))+Object$Pred_T 
    CIs_UB <- (crit_val*sqrt(Object$Var_Delta.T_S)) +Object$Pred_T 
    
    aantal_Zero_included <- length(Object$Var_Delta.T[CIs_LB<0 & CIs_UB>0])
    aantal_CI_under_Zero <- length(Object$Var_Delta.T[CIs_LB<0 & CIs_UB<0])
    aantal_CI_boven_Zero <- length(Object$Var_Delta.T[CIs_LB>0 & CIs_UB>0])
    
    temp_vect <- rep(NA, length(Object$Var_Delta.T_S))
    temp_vect[CIs_LB<0 & CIs_UB>0] <- 1 #zero included
    temp_vect[CIs_LB<0 & CIs_UB<0] <- 2  # under zero
    temp_vect[CIs_LB>0 & CIs_UB>0] <- 3  # above zero
    
    if (aantal_Zero_included>0){
      cat("0 included in support interval: ", aantal_Zero_included/length(Object$Var_Delta.T_S),
          "     (obtained for rho_T0T1 values in range [", min(Object$T0T1[temp_vect==1]), "; ", 
          max(Object$T0T1[temp_vect==1]), "])\n", sep="")
    }
    if (aantal_Zero_included==0){
      cat("0 included in support interval: ", aantal_Zero_included/length(Object$Var_Delta.T_S), "\n")
    }  
    
    if (aantal_CI_under_Zero>0){
      cat("Entire support interval below 0: ", aantal_CI_under_Zero/length(Object$Var_Delta.T_S),
          "     (obtained for rho_T0T1 values in range [", min(Object$T0T1[temp_vect==2]), "; ", 
          max(Object$T0T1[temp_vect==2]), "])\n", sep="")
    }
    if (aantal_CI_under_Zero==0){
      cat("Entire support interval below 0: ", aantal_CI_under_Zero/length(Object$Var_Delta.T_S), "\n")
    }
    
    if (aantal_CI_boven_Zero > 0){
      cat("Entire support interval above 0: ", aantal_CI_boven_Zero/length(Object$Var_Delta.T_S),
          "     (obtained for rho_T0T1 values in range [", min(Object$T0T1[temp_vect==3]), "; ", 
          max(Object$T0T1[temp_vect==3]), "])\n", sep="")
    }
    if (aantal_CI_boven_Zero == 0){
      cat("Entire support interval above 0: ", aantal_CI_boven_Zero/length(Object$Var_Delta.T_S), "\n")
    }
  }
  
  

}