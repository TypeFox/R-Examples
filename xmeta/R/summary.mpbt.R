summary.mpbt <- function(object,...) {

summary_mpbt <- function(object) {
  type <- object$type
  k <- object$k
  method <- object$method

  cat("Outcome:",type,fill=TRUE)
  if (method=="nn.cl") {
    cat("Method: score test for detecting publication bias when the within-study correlations are unknown",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }

 

  ## print the object
 
  if(method=="nn.cl"){
  cat("Test statistics (mpbt):",object$mpbt.TS, fill=TRUE)
  

  cat("P value:",object$mpbt.pv, fill=TRUE)
 
  cat("\n")
}

}

  if (!inherits(object, "mpbt"))
    stop("Use only with 'mpbt' objects.\n")
	result <- summary_mpbt(object)
  invisible(result)
}



