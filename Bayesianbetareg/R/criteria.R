criteria <-
function(X,beta.residuals){
  if(is.null(beta.residuals)){
    stop("Residual data is not included")
  }
  if(is.null(X)){
    stop("Variable data is not included")
  }
  deviance = sum(beta.residuals$deviance^2)
  AIC <- deviance + (2*(dim(X)[2]))  
  BIC <- deviance+(log(dim(X)[1])*(dim(X)[2]))
  
  criteria <- list()
  
  criteria$Deviance <- deviance
  criteria$AIC <- AIC
  criteria$BIC <- BIC
  
  criteria
}
