##' @export
print.splitMethod <- function(x,...){
  if (x$name=="no plan")
    return(cat("\nNo data splitting: either apparent or independent test sample performance\n"))
  cat("\nMethod for estimating the prediction error:\n")
  if (x$internal.name=="crossval"){
    cat("\n",x$name,"\n\n")
    cat("Repeated: ",x$B,ifelse(x$B==1," time","times"),"\n")
  }
  else{
    if (x$internal.name=="loocv"){
      cat("\n",x$name,"\n\n")
    }
    else{
      cat("\nBootstrap cross-validation\n\n")
      if (x$M<x$N)
        cat("Type: subsampling\n",x$M," out of ",x$N,"\n\n")
      else
        cat("Type: resampling\nBootstrap sample size: ",x$N,"\n")
      cat("No. bootstrap samples: ",x$B,"\n")
    }
    cat("Sample size: ",x$N,"\n")
  }
}
