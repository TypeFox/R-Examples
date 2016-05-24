summary.galaxy <- function(object,...) {

summary_galaxy <- function(object) {
  type <- object$type
  k <- object$k
  method <- object$method
  L <- object$L

  cat("Outcome:",type,fill=TRUE)
  if (method=="galaxy.cl") {
    cat("Method: galaxy method for adjusting publication bias when the within-study correlations are unknown",fill=TRUE)
    cat("Number of outcomes:",k,"\n")
  }

 

  ## print the object
 
  if(method=="galaxy.cl"){
  cat("Coefficients:",object$res_galaxy, fill=TRUE)  
  cat("Variances:",object$cov_galaxy, fill=TRUE)
  cat("k0:",object$k0, fill=TRUE)
  cat("side:",object$side, fill=TRUE)
  cat("\n")
}

}

  if (!inherits(object, "galaxy"))
    stop("Use only with 'galaxy' objects.\n")
	result <- summary_galaxy(object)
  invisible(result)
}



