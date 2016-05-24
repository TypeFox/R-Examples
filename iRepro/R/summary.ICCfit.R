summary.ICCfit <-
function(object, ...){
  t <- as.data.frame(c(object$icc,object$sigma2.b,object$sigma2.w,object$mu))
  rownames(t) <- c("ICC", "Between-class variance", "Within-class variance","Mu")
  colnames(t) <- "Estimate"
  
  loglik <- as.data.frame(object$loglikelihood)
  colnames(loglik) <- ""
  rownames(loglik) <- "Log-likelihood: "
  
  res <- list(estimates = t,
              loglikelihood = loglik)
  
  class(res) <- "summary.ICCfit"
  res
}

