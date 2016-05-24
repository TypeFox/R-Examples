### summary function for geem object.
summary.geem <- function(object, ...)  {
  Coefs <- matrix(NA,nrow=length(object$beta),ncol=5)
  Coefs[,1] <- c(object$beta)
  naive <- is.character(object$var)
  if(!naive & any(diag(object$var) < 0) ){
    naive <- TRUE
    warning("Some elements of robust variance estimate < 0.  Reporting model based SE.")
  }
  Coefs[,2] <- sqrt(diag(object$naiv.var))
  if(naive){Coefs[,3] <- rep(NA, length(object$beta))}else{Coefs[,3] <- sqrt(diag(object$var))}
  if(naive){Coefs[,4] <- Coefs[,1]/Coefs[,2]}else{Coefs[,4] <- Coefs[,1]/Coefs[,3]}
  Coefs[,5] <- round(2*pnorm(abs(Coefs[,4]), lower.tail=F), digits=8)
  colnames(Coefs) <- c("Estimates","Model SE","Robust SE", "wald", "p")
  
  summ <- list(beta = Coefs[,1], se.model = Coefs[,2], se.robust = Coefs[,3], wald.test = Coefs[,4], p = Coefs[,5],
               alpha = object$alpha, corr = object$corr, phi = object$phi, niter = object$niter, clusz = object$clusz,
               coefnames = object$coefnames, weights=object$weights, biggest.R.alpha = object$biggest.R.alpha)
  
  class(summ) <- 'summary.geem'
  return(summ)
}



### print function for summary.geem object
print.summary.geem <- function(x, ...){
  Coefs <- matrix(0,nrow=length(x$coefnames),ncol=5)
  rownames(Coefs) <- c(x$coefnames)
  colnames(Coefs) <- c("Estimates","Model SE","Robust SE", "wald", "p")
  Coefs[,1] <- x$beta
  Coefs[,2] <- x$se.model
  Coefs[,3] <- x$se.robust
  Coefs[,4] <- x$wald.test
  Coefs[,5] <- x$p
  
  #print("Call: ", object$call, "\n")
  print(signif(Coefs, digits=4))
  if(!is.element(x$corr, c("independence", "ar1", "exchangeable"))){
    if(dim(x$biggest.R.alpha)[1] > 4 ){
      cat("\n Working Correlation[1:4,1:4]: \n")
      print(as.matrix(round(x$biggest.R.alpha[1:4,1:4], digits=4)))
    }else{
      cat("\n Working Correlation: \n")
      print(as.matrix(round(x$biggest.R.alpha, digits=4)))
    }
  }else{
    cat("\n Estimated Correlation Parameter: ", signif(x$alpha, digits=4), "\n")
  }
  #cat("\n Est. Correlation: ", signif(x$alpha, digits=4), "\n")
  cat(" Correlation Structure: ", x$corr, "\n")
  cat(" Est. Scale Parameter: ", signif(x$phi, digits=4), "\n")
  cat("\n Number of GEE iterations:", x$niter, "\n")
  cat(" Number of Clusters: ", length(x$clusz), "   Maximum Cluster Size: ", max(x$clusz), "\n")
  cat(" Number of observations with nonzero weight: ", sum(x$weights != 0), "\n")
}

### print function for geem object
print.geem <- function(x, ...){
  coefdf <- signif(data.frame(x$beta), digits=4)
  rownames(coefdf) <- x$coefnames
  colnames(coefdf) <- ""
  print(x$call)
  cat("\n", "Coefficients:", "\n")
  print(t(coefdf))
  cat("\n Scale Parameter: ", signif(x$phi, digits=4), "\n")
  cat("\n Correlation Model: ", x$corr)
  if(!is.element(x$corr, c("independence", "ar1", "exchangeable"))){
    if(dim(x$biggest.R.alpha)[1] > 4 ){
      cat("\n Working Correlation[1:4,1:4]: \n")
      print(as.matrix(round(x$biggest.R.alpha[1:4,1:4], digits=4)))
    }else{
      cat("\n Working Correlation: \n")
      print(as.matrix(round(x$biggest.R.alpha, digits=4)))
    }
  }else{
    cat("\n Estimated Correlation Parameter: ", signif(x$alpha, digits=4), "\n")
  }
  cat("\n Number of clusters: ", length(x$clusz), "  Maximum cluster size: ", max(x$clusz), "\n")
  cat(" Number of observations with nonzero weight: ", sum(x$weights != 0), "\n")
}

#print.coef.geem <- function(x, ...){
#  print(signif(x, 3))
#}