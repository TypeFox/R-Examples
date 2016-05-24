summary.mmeta <- function(object,...) {

summary_mmeta <- function(object) {
  type <- object$type
  k <- object$k
  method <- object$method

  cat("Outcome:",type,fill=TRUE)
  if (method=="nn.reml") {
    cat("Method: restricted maximum likelihood method",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }

  if (method=="nn.cl") {
    cat("Method: composite likelihood method",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }

   if (method=="nn.mom") {
    cat("Method: method of moment",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }

   if (method=="nn.mom") {
    cat("Method: improved method for Riley model",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }

   if (method=="bn.cl") {
    cat("Method: marginal bivariate normal model",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }


   if (method=="bb.cl") {
    cat("Method: marginal beta-binomial model",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }
   if (method=="tb.cl") {
    cat("Method: hybrid model for disease prevalence along with sensitivity and specificity for diagnostic test accuracy",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }

   if (method=="tn.cl") {
    cat("Method: trivariate model for multivariate meta-analysis of diagnostic test accuracy",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }

  ## print the object
  if(method=="nn.reml"){
  cat("Fixed-effects coefficients:",fill=TRUE)
  beta.coefficients=data.frame(cbind(object$coefficients[1],object$coefficients[3]))
  colnames(beta.coefficients)=c("beta1", "beta2")
  rownames(beta.coefficients)=c("")
  print(beta.coefficients)

  cat("Between-study variances:",fill=TRUE)
  tau2.coefficients=data.frame(cbind(object$coefficients[2],object$coefficients[4]))
  colnames(tau2.coefficients)=c("tau1.2", "tau2.2")
  rownames(tau2.coefficients)=c("")
  print(tau2.coefficients)

  cat("Between-study correlation:",fill=TRUE)
  rhob=data.frame(object$coefficients[5])
  colnames(rhob)=c("rhob")
  rownames(rhob)=c("")
  print(rhob)
  cat("Covariance matrix:",fill=TRUE)
  print(object$vcov)
  cat("\n")
}

  if(method=="nn.cl"){
  cat("Fixed-effects coefficients:",fill=TRUE)
  beta.coefficients=data.frame(cbind(object$coefficients[1],object$coefficients[3]))
  colnames(beta.coefficients)=c("beta1", "beta2")
  rownames(beta.coefficients)=c("")
  print(beta.coefficients)

  cat("Between-study variances:",fill=TRUE)
  tau2.coefficients=data.frame(cbind(object$coefficients[2],object$coefficients[4]))
  colnames(tau2.coefficients)=c("tau1.2", "tau2.2")
  rownames(tau2.coefficients)=c("")
  print(tau2.coefficients)
  cat("Covariance matrix:",fill=TRUE)
  print(object$vcov)
  cat("\n")
}

  if(method=="nn.mom"){
  cat("Fixed-effects coefficients:",fill=TRUE)
  beta.coefficients=data.frame(cbind(object$coefficients[1],object$coefficients[2]))
  colnames(beta.coefficients)=c("beta1", "beta2")
  rownames(beta.coefficients)=c("")
  print(beta.coefficients)
  cat("Covariance matrix:",fill=TRUE)
  print(object$vcov)
  cat("\n")
}

  if(method=="nn.rs"){
   cat("Fixed-effects coefficients:",fill=TRUE)
  beta.coefficients=data.frame(cbind(object$coefficients[1],object$coefficients[3]))
  colnames(beta.coefficients)=c("beta1", "beta2")
  rownames(beta.coefficients)=c("")
  print(beta.coefficients)

  cat("Between-study variances:",fill=TRUE)
  tau2.coefficients=data.frame(cbind(object$coefficients[2],object$coefficients[4]))
  colnames(tau2.coefficients)=c("tau1.2", "tau2.2")
  rownames(tau2.coefficients)=c("")
  print(tau2.coefficients)

  cat("Between-study correlation:",fill=TRUE)
  rhos=data.frame(object$coefficients[5])
  colnames(rhos)=c("rhos")
  rownames(rhos)=c("")
  print(rhos)
  cat("Covariance matrix:",fill=TRUE)
  print(object$vcov)
  cat("\n")
}

if(method=="bb.cl"){
   cat("Log odds ratio:", object$logOR, fill=TRUE)
   cat("Ccoffidence interval for odds ratio:",object$OR_CI,fill=TRUE)
   cat("Standard error for log odds ratio:",object$se,fill=TRUE) 
   cat("Hessian matrix in log scale:",fill=TRUE) 
   print(object$hessian.log)
   cat("Convergence:",object$conv, fill=TRUE)
   cat("MLE",fill=TRUE) 
   mle=data.frame(t(object$mle[c(1:4)]))
   colnames(mle)=c("a1", "b1","a2", "b2")
   rownames(mle)=c("")
   print(mle)
   cat("\n")
}

 if(method=="bn.cl"){
  cat("Coefficients:", fill=TRUE)
  print(object$coefficients)
  cat("Covariance matrix:",fill=TRUE)
  print(object$vcov)
  cat("\n")
}
 if(method=="tb.cl"){
  cat("Coefficients:", fill=TRUE)
  coefficients=t(object$coefficients)
  rownames(coefficients)=""
  print(coefficients)
  cat("Covariance matrix:",fill=TRUE)
  print(object$vcov)
  cat("\n")
}

 if(method=="tn.cl"){
  cat("Coefficients:", fill=TRUE)
  coefficients=t(object$coefficients)
  rownames(coefficients)=""
  print(coefficients)
  cat("Covariance matrix:",fill=TRUE)
  print(object$vcov)
  cat("\n")
}
}

  if (!inherits(object, "mmeta"))
    stop("Use only with 'mmeta' objects.\n")
	result <- summary_mmeta(object)
  invisible(result)
}



