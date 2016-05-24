summary.propagate <- function(object, ...)
{
  ## print error propagation results
  cat("Results from error propagation:\n")
  print(object$prop)
  
  ## print simulation results
  cat("\nResults from Monte Carlo simulation:\n")
  print(object$sim)
  
  ## print covariance matrix
  cat("\nCovariance matrix:\n")
  print(object$covMat)
  
  ## print derivatives of gradient
  cat("\nSymbolic gradient matrix:\n")
  print(as.character(object$gradient))
  
  ## print derivatives of gradient
  cat("\nEvaluated gradient matrix:\n")
  print(object$evalGrad)
   
  ## print derivatives of hessian
  cat("\nSymbolic hessian matrix:\n")
  NCOL <- ncol(object$hessian)    
  STR <- as.character(object$hessian) 
  SPLIT <- split(STR, as.factor(rep(1:NCOL, NCOL)))
  for (i in 1:length(SPLIT)) print(SPLIT[[i]])
    
  ## print derivatives of hessian
  cat("\nEvaluated hessian matrix:\n")
  NCOL <- ncol(object$hessian)
  EVAL <- object$evalHess 
  SPLIT <- split(EVAL, as.factor(rep(1:NCOL, NCOL)))
  for (i in 1:length(SPLIT)) print(SPLIT[[i]])  
  
  ## Skewness and excess kurtosis of evaluated MC simulations
  if (length(object$resSIM) > 1) {
    cat("\nSkewness / Excess Kurtosis of MC evaluations:\n")
    cat(skewness(object$resSIM), "/", kurtosis(object$resSIM), "\n")
  }
    
  ## Shapiro test for normality of MC distribution
  if (length(object$resSIM) > 1) {
    cat("\nShapiro-Wilk test for normality: ")
    if (length(object$resSIM) > 5000) 
    DAT <- object$resSIM[1:5000] else DAT <- object$resSIM 
    PVAL <- shapiro.test(DAT)$p.value
    cat(PVAL)
    if (PVAL >= 0.05) cat(" => normal") else cat(" => non-normal")
  }
    
  ## Kolmogorov-Smirnov test for normality of MC distribution
  if (length(object$resSIM) > 1) {
    cat("\nKolmogorov-Smirnov test for normality: ")
    DAT <- object$resSIM
    simDAT <- rnorm(length(DAT), mean(DAT, na.rm = TRUE), sd(DAT, na.rm = TRUE))
    PVAL <- ks.test(DAT, simDAT)$p.value
    cat(PVAL)
    if (PVAL >= 0.05) cat(" => normal") else cat(" => non-normal")  
  }
}
