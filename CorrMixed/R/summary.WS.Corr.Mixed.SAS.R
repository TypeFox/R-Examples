summary.WS.Corr.Mixed.SAS <- function(object, ..., Object){

  if (missing(Object)){Object <- object} 
  x <- Object
  
  if (Object$Model=="Model 1"){
    cat(Object$Model, ", Random intercept:\n", sep="")    
    cat("==========================")
    cat("\n\nEstimated correlations R (r(time_j, time_k) constant): \n")
    cat("------------------------------------------------------ \n")  
    cat("R:", Object$R[1])
    cat("\n\n", (1-Object$Alpha)*100, "% confidence interval (Delta method): [", 
        Object$CI.Lower[1], "; ", Object$CI.Upper[1], "]", sep="")
  }
  
  if (Object$Model=="Model 2"){
    cat(Object$Model, ", Random intercept + serial corr (Gaussian):\n", sep="")    
    cat("===================================================\n\n")
    cat("Estimated correlations R as a function of time lag: \n")
    cat("--------------------------------------------------- \n")  
    print(Object$R)
    cat("\n", (1-Object$Alpha)*100, "% confidence intervals (Delta method), lower bounds:\n", sep="")
    cat("------------------------------------------------------ \n")
    print(Object$CI.Lower)
    cat("\n", (1-Object$Alpha)*100, "% confidence intervals (Delta method), upper bounds:\n", sep="")
    cat("------------------------------------------------------ \n")
    print(Object$CI.Upper)
  }
  
  if (Object$Model=="Model 3"){
    cat(Object$Model, ", Random intercept, slope + serial corr (Gaussian):\n", sep="")    
    cat("==========================================================\n\n")
    cat("\n\nEstimated correlations R at each time point r(time_j, time_k) \n")
    cat("------------------------------------------------------------- \n")  
    print(Object$R)
    cat("\n\nEstimated correlations R as a function of time lag (loess) \n")
    cat("---------------------------------------------------------- \n")  
    print(Object$Pred.Model3.Loess)
  }
}
