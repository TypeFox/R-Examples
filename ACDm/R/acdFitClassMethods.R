coef.acdFit <- function(object, returnCoef = "all", ...){
  returnCoef <- match.arg(returnCoef, c("all", "distribution", "model"))
  
  switch(returnCoef,
                        all = c(object$mPara, object$dPara),
                        distribution = object$dPara,
                        model = object$mPara)
}

residuals.acdFit <- function(object, ...){
  object$residuals
}

print.acdFit <- function(x, ...){
  
  if(x$distribution == "exponential") {
    cat("\nACD model estimation by (Quasi) Maximum Likelihood \n")
  } else {
    cat("\nACD model estimation by Maximum Likelihood \n")
  }
  
  cat("\nCall:\n") 
  cat(" ", deparse(x$call), "\n")
  cat("\nModel:\n") 
  cat(" ", x$model)
  cat("(")
  cat(x$order[1])
  for(i in 2:length(x$order)) cat("", x$order[i], sep = ", ")
  cat(")")
  if(length(x$SNIACDbp) != 0) cat("\n  Break points:", x$SNIACDbp)
  cat("\n")
  cat("\nDistribution:\n") 
  cat(" ", x$distribution)
  cat("\n\nN:", x$N) 
  cat("\n\nParameter estimate:\n") 
  print(format(x$parameterInference, digits = 3, scientific = F))
  if(length(x$comments) > 0){
    cat("\nNote:", x$comments) 
  }
  if(length(x$forcedDistPara) > 0){
    cat("\n\nThe fixed/unfree mean distribution parameter: \n") 
    cat(" ", names(x$forcedDistPara), ": ", x$forcedDistPara, sep = "")
  }
  if(length(x$bootErr) != 0){
    cat("\n\nBootstrap correlations:\n")
    print(format(data.frame(x$bootCorr), digits = 3, scientific = F))
  }
  if(length(x$robustCorr) != 0){
    cat("\n\nQML robust correlations:\n")
    print(format(data.frame(x$robustCorr), digits = 3, scientific = F))
  }
  cat("\n\nGoodness of fit:\n")
  print.data.frame(x$goodnessOfFit)
  cat("\nConvergence:", x$convergence, "\n")  
  cat("\nNumber of log-likelihood function evaluations:", x$evals, "\n")  
  if(length(x$bootErr) == 0) cat("\nEstimation time:", round(x$estimationTime, digits = 4), attributes(x$estimationTime)$units, "\n") 
  else cat("\nTotal estimation time (including bootstrap simulations):", round(x$estimationTime, digits = 4), attributes(x$estimationTime)$units, "\n") 
  cat("\nDescription:", x$description)
  cat("\n\n")
}

