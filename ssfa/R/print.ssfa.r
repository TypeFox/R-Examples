print.ssfa <- function(x, ...) {
  
  if(x$rho==0)
  {    
  cat("\n")
  cat("Stochastic frontier analysis model\n")
  cat("\n")
  cat("Coefficients:\n")
    
  coef <-c(x$coef["Intercept"], 
           x$coef[names(x$coef)[4:length(x$coef)]], 
           x$coef["sigmau2"], 
           x$coef["sigmav2"], 
           x$sigma2)
  names(coef) <- c("Intercept", 
                          names(x$coef)[4:length(x$coef)], 
                          "sigmau2", 
                          "sigmav2", 
                          "sigma2")
  }
  
  if(x$rho!=0)
  {    
    cat("\n")
    cat("Spatial Stochastic frontier analysis model\n")
    cat("\n")
    cat("Pay attention:")
    cat("\n")
    cat("1 - classical SFA sigmau2 = sigmau2_dmu + sigmau2_sar")          
    cat("\n")
    cat("2 - sigma2 = sigmau2_dmu + sigmau2_sar + sigmav2")
    cat("\n")
    cat("\n")
    cat("Coefficients:\n")
    cat("\n")
    coef <- c(x$coef["Intercept"], 
              x$coef[names(x$coef)[5:length(x$coef)]],
              x$coef["sigmau2_dmu"],
              x$coef["sigmav2"], 
              x$sigmau2_sar, 
              x$sigma2, 
              x$coef["rho"])
    names(coef) <- c("Intercept", 
                            names(x$coef)[5:length(x$coef)], 
                            "sigmau2_dmu", 
                            "sigmav2", 
                            "sigmau2_sar", 
                            "sigma2", 
                            "rho")
  }
  
  print(coef)
  cat("\n")
  invisible(x)
}