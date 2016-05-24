`summary.ppar` <-
function(object,...)
# summary method for object of class "ppar"
{
  
  if (length(object$pers.ex) > 0) {
    thetaind <- rownames(object$X)[-object$pers.ex]
  } else {
    thetaind <- rownames(object$X)
  }
    
  if (any(is.na(object$X))) {                                       #recompute gmemb without persons excluded
    dichX <- ifelse(is.na(object$X),1,0)
    strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
    gmemb <- as.vector(data.matrix(data.frame(strdata)))
  } else {
    gmemb <- rep(1,dim(object$X)[1])
  }
  
  cat("\n")
  cat("Estimation of Ability Parameters")
  for (i in 1:length(object$thetapar)) {
    cat("\n\n")
    if (length(object$thetapar) > 1) {
      cat("Subject NA Group:",i,"\n")
      xvec <- rbind(object$X[gmemb==i,])[1,]                    #determine NA pattern
      xvec[!is.na(xvec)] <- "x"
      cat("NA pattern:",xvec,"\n")
      }
    cat("Collapsed log-likelihood:",object$loglik[[i]],"\n")
    cat("Number of iterations:",object$iter[[i]],"\n")
    cat("Number of parameters:",object$npar[[i]],"\n")
    cat("\n")
    cat("ML estimated ability parameters (without spline interpolated values): \n")
    coef.table <- cbind(object$thetapar[[i]],object$se.theta[[i]],confint(object)[[i]])
    dimnames(coef.table) <- list(paste("theta",thetaind[object$gmemb==i]),c("Estimate","Std. Err.",colnames(confint(object)[[i]])))
    print(coef.table)
  }
}

