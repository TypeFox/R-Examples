
summary.mcmc_hsar <- function(object, ...)
{
  x<-object
  cat("\nCall:\n")
  if ( class(x)=="mcmc_hsar" ){
    cat("hsar( X, y, W, M, Z, Unum ) \n")
  }
  if ( class(x)=="mcmc_hsar_rho_0" ){
    cat("hsar( X, y, M, Z, Unum ) \n")
  }
  if ( class(x)=="mcmc_hsar_lambda_0" ){
    cat("hsar( X, y, W, Z, Unum ) \n")
  }
  
  cat("\nCoefficients:\n")
  print(x$Mbetas)
  cat("\n")
  
  cat("\n Spatial coefficients\n")
  cat("Strength of the spatial interaction rho",x$Mrho, "\n")
  cat("Strength of the spatial interaction at higher level lambda",x$Mlambda, "\n")
  
  cat("\n Diagnostics \n")
  cat("\n Deviance information criterion (DIC):", x$DIC, "\n")
  cat("Effective number of parameters (pd):", x$pD, "\n")
  cat("\n Log likelihood:", x$Log_Likelihood, "\n")
  cat("\n Pseudo R squared:", x$R_Squared, "\n")
  
  if ( class(x)=="mcmc_hsar" ){
    cat("\n Total impact: \n",x$impact_total, "\n")
    cat("Indirect impact: \n",x$impact_idirect, "\n")
    cat("Direct impact: \n",x$impact_direct, "\n")
  }
  invisible(x)
}

summary.mcmc_sar <- function(object, ...)
{
  x<-object
  cat("\nCall:\n")
  cat("sar( X, y, W ) \n")
  
  cat("\nCoefficients:\n")
  print(x$Mbetas)
  cat("\n")
  
  cat("\n Spatial coefficients\n")
  cat("Strength of the spatial interaction rho",x$Mrho, "\n")
  
  cat("\n Diagnostics \n")
  cat("\n Deviance information criterion (DIC):", x$DIC, "\n")
  cat("Effective number of parameters (pd):", x$pD, "\n")
  cat("\n Log likelihood:", x$Log_Likelihood, "\n")
  cat("\n Pseudo R squared:", x$R_Squared, "\n")
  cat("\n Total impact: \n",x$impact_total, "\n")
  cat("Indirect impact: \n",x$impact_idirect, "\n")
  cat("Direct impact: \n",x$impact_direct, "\n")
  invisible(x)
}

summary.mcmc_hsar_rho_0 <- function(object, ...)
{
  x <- object 
  cat("\nCall:\n")
  cat("hsar( X, y, M, Z, Unum ) \n")
  
  cat("\nCoefficients:\n")
  print(x$Mbetas)
  cat("\n")
  
  cat("\n Spatial coefficients\n")
  cat("Strength of the spatial interaction at higher level lambda",x$Mlambda, "\n")
  
  cat("\n Diagnostics \n")
  cat("\n Deviance information criterion (DIC):", x$DIC, "\n")
  cat("Effective number of parameters (pd):", x$pD, "\n")
  cat("\n Log likelihood:", x$Log_Likelihood, "\n")
  cat("\n Pseudo R squared:", x$R_Squared, "\n")
  invisible(x)
}

summary.mcmc_hsar_lambda_0 <- function(object, ...)
{
  x <- object
  cat("\nCall:\n")
  cat("hsar( X, y, W, Z, Unum ) \n")
  
  cat("\nCoefficients:\n")
  print(x$Mbetas)
  cat("\n")
  
  cat("\n Spatial coefficients\n")
  cat("Strength of the spatial interaction rho",x$Mrho, "\n")
  
  cat("\n Diagnostics \n")
  cat("\n Deviance information criterion (DIC):", x$DIC, "\n")
  cat("Effective number of parameters (pd):", x$pD, "\n")
  cat("\n Log likelihood:", x$Log_Likelihood, "\n")
  cat("\n Pseudo R squared:", x$R_Squared, "\n")
  invisible(x)
}

