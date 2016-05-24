print.cv.BTLLasso <- function(x, ...){
  
  m <- x$m
  n.theta <- x$n.theta
  
  covar <- colnames(x$X)
  
  cat("Output of BTL-Lasso estimation:","\n")
  
  cat("---","\n")

  cat("Setting:")
  cat("\n", x$n, "persons")
  cat("\n", x$m, "items")
  cat("\n", x$q+1, "response categories")
  cat("\n", x$p, "covariates")
  cat("\n", length(x$lambda), "different tuning parameters","\n")
  
  cat("---","\n")
  
  cat("Parameter estimates after",x$folds,"-","fold cross-validation","\n")
  
  coefs <- x$coefs.repar[which.min(x$deviance),]
  theta <- coefs[1:n.theta]
  
  gamma0 <- coefs[(n.theta+1):(n.theta+m)]
  
  gamma.mat <- matrix(coefs[(n.theta+m+1):length(coefs)], byrow=TRUE, ncol=m)
  
  coef.mat <- rbind(gamma0,gamma.mat)
  
  rownames(coef.mat) <- c("(Intercept)",covar)
  colnames(coef.mat) <- x$labels
  

  print(coef.mat, ...)
  
  cat("\n")
  
  cat("Optimal lambda:",x$lambda[which.min(x$deviance)],"\n")
  
  cat("\n")
  
  cat("log likelihood:",x$logLik[which.min(x$deviance)],"\n")
  
  invisible(x)
  
  
}

