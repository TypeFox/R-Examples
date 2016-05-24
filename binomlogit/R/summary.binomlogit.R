summary.binomlogit <- function(object, ...){
  cat("\nData input:\n\n")
  if(is(object,"binomlogitIndiv")){
    cat("Observations:", object$N)
  } else {
    cat("Observations:", object$t)
  }
  cat("\nCovariates:", object$dims)
  
  cat("\n\n\nMCMC details:\n\n")
  cat("Draws (incl. BI):", object$sim, "draws\n")
  cat("Burn-in:", object$burn, "draws\n")
  if(is(object,"binomlogitMH")){
    cat("Accept. phase during BI:", object$acc, "draws\n")
  }
  if(is(object,"binomlogitHAM")){
    cat("Boundaries for HAM:", object$low,"/",object$up,"\n")
  }
  
  cat("\n\nPrior parameters:\n\n")
  cat("b0\n")
  print(object$b0)
  cat("\nB0\n")
  print(object$B0)
  
  cat("\n\nRuntime:\n\n")
  cat("Total time: ", object$duration, " sec.\n")
  cat("Time (without BI):", object$duration_wBI, "sec.\n")
  
  if(is(object,"binomlogitMH")||is(object,"binomlogitHAM")){
    cat("\n\nAcceptance rate:", object$rate, "%\n")
  }
  
  mw_beta=as.matrix(sapply(1:object$dims,function(y) mean(object$beta[y,(object$burn+1):object$sim])))
  cat("\n\nPosterior mean:\n\n")
  cat("beta\n")
  print(mw_beta)
}
