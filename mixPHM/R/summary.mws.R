`summary.mws` <-
function(object,...)
#summary method for object of class "mws"
{
  cat("Likelihood value:",object$likelihood[length(object$likelihood)],"\n")
  cat("Number of estimated parameters:",object$npar,"\n")
  cat("AIC:",object$bic,"\n")
  cat("BIC:",object$aic,"\n")
  cat("\n")
  cat("\n")
  cat("Scale parameters: \n")
  rownames(object$scale) <- paste("Component",1:object$K,sep="")
  colnames(object$scale) <- paste("V",1:(dim(object$scale)[2]),sep="")
  print(object$scale)
  cat("\n")
  cat("Shape parameters: \n")
  rownames(object$shape) <- paste("Component",1:object$K,sep="")
  colnames(object$shape) <- paste("V",1:(dim(object$shape)[2]),sep="")
  print(object$shape)
  cat("\n")
}

