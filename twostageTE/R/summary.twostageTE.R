#' Summary twostageTE object.
#' 
#' Summary of twostageTE object
#' 
## #' @usage summary(object, ...)
## #' @aliases summary summary.twostageTE
#' @param object twostageTE object
#' @param ... ignored
#' @S3method print twostageTE
summary.twostageTE <-
function(object,...){	
	if(!inherits(object,"twostageTE")){
		stop("Error:  Object is not of class twostageTE")
	}

	if(!is.null(cl <- object$call)) {
	names(cl)[2] <- ""
	cat("Call:\n")
	dput(cl)
	}

	cat(sprintf("\n%.1f%% Confidence Interval", object$level*100))
	if (is.na(object$L2)) {
		cat(sprintf("\nn   Lower   d0_hat   Upper\n%d   %.3f   %.3f   %.3f\n",length(object$Y1),object$L1,object$estimate,object$U1))
		if (length(grep(pattern="IR-likelihood", object$call)) == 0) {
			cat(sprintf("\nAuxiliary Estimates\nf'(d_0)   sigma^2\n%.3f     %.3f\n",object$deriv_d0,object$sigmaSq))
		}
		else {
			cat(sprintf("\nAuxiliary Estimates\nf'(d_0)   sigmaSq\nNA     %.3f\n",object$sigmaSq))
		}
	}
	else {
		cat(sprintf("\nn1   n2   Lower   d0_hat   Upper\n%d   %d   %.3f   %.3f   %.3f\n",length(object$Y1),length(object$Y2),object$L2,object$estimate,object$U2))
		if (length(grep(pattern="IR-likelihood", object$call)) == 0) {
			cat(sprintf("\nAuxiliary Estimates\nf'(d_0)   sigma^2\n%.3f     %.3f\n",object$deriv_d0,object$sigmaSq))
		}
		else {
			cat(sprintf("\nAuxiliary Estimates\nf'(d_0)   sigmaSq\nNA     %.3f\n",object$sigmaSq))
		}
	}
  	invisible(object)
}
