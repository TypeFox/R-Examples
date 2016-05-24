summary.WS.Corr.Mixed <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  x <- Object
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n")
  
  if (Object$Model=="Model 1, Random intercept"){
    cat(Object$Model, "\n")    
    cat("=========================\n\n")
    cat("Fitted variance components: \n")
    cat("--------------------------- \n")
    cat("D:", Object$D, "\n")
    cat("Sigma**2:", Object$Sigma2, "\n")
    cat("\nEstimated correlations R (r(time_j, time_k) constant): \n")
    cat("------------------------------------------------------ \n")  
    cat("R: ", (Object$R)[1], sep="")
    cat("\n", (1-Object$Alpha)*100, "% confidence interval (bootstrap): [", 
        Object$CI.Lower[1], "; ", Object$CI.Upper[1], "]", sep="")
    cat("\n\nModel fit: \n")
    cat("---------- \n")  
    cat("LogLik: ", Object$LogLik)
    cat("\nAIC: ", Object$AIC)
  }
  
  if (Object$Model=="Model 2, Random intercept + serial corr (Gaussian)"){
    cat(Object$Model, "\n")    
    cat("==================================================\n\n")
    cat("Fitted variance components: \n")
    cat("--------------------------- \n")
    cat("D:", Object$D, "\n")
    cat("Sigma**2:", Object$Sigma2, "\n")
    cat("Tau**2:", Object$Tau2, "\n")
    cat("Rho:", Object$Rho, "\n")  
    cat("\nEstimated correlations R as a function of time lag: \n")
    cat("--------------------------------------------------- \n")  
    print(Object$R)
    cat("\n", (1-Object$Alpha)*100, "% confidence intervals (bootstrap), lower bounds:\n", sep="")
    cat("--------------------------------------------------- \n")  
    print(Object$CI.Lower)
    cat("\n", (1-Object$Alpha)*100, "% confidence intervals (bootstrap), upper bounds:\n", sep="")
    cat("--------------------------------------------------- \n")  
    print(Object$CI.Upper)
    cat("\n\nModel fit: \n")
    cat("---------- \n")  
    cat("LogLik: ", Object$LogLik)
    cat("\nAIC: ", Object$AIC)
  }
  if (Object$Model=="Model 3, Random intercept, slope + serial corr (Gaussian)"){
    cat(Object$Model, "\n")    
    cat("=========================================================\n\n")
    cat("Fitted variance components: \n")
    cat("--------------------------- \n")
    cat("D:\n")
    print(Object$D)
    cat("\nSigma**2:", Object$Sigma2, "\n")
    cat("Tau**2:", Object$Tau2, "\n")
    cat("Rho:", Object$Rho, "\n")  
    cat("\nEstimated correlations R at each time point r(time_j, time_k) \n")
    cat("------------------------------------------------------------- \n")  
    print(Object$R)
    cat("\n", (1-Object$Alpha)*100, "% confidence intervals (bootstrap), lower bounds:\n", sep="")
    cat("--------------------------------------------------- \n") 
    print(Object$CI.Lower)
    cat("\n", (1-Object$Alpha)*100, "% confidence intervals (bootstrap), upper bounds:\n", sep="")
    cat("--------------------------------------------------- \n") 
    print(Object$CI.Upper)
    cat("\n\nModel fit: \n")
    cat("---------- \n")  
    cat("LogLik: ", Object$LogLik)
    cat("\nAIC: ", Object$AIC)
  }
}
