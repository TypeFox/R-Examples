#' Prints twostageTE object.
#' 
#' Prints twostageTE object
#' 
## #' @usage print(x, ...)
## #' @aliases print print.twostageTE
#' @param x twostageTE object
#' @param ... ignored
#' @S3method print twostageTE
print.twostageTE <- function(x, ...){	
	if(!inherits(x,"twostageTE")){
		stop("Error:  Object is not of class twostageTE")
	}

	if(!is.null(cl <- x$call)) {
	names(cl)[2] <- ""
	cat("Call:\n")
	dput(cl)
	}

	cat(sprintf("\n%.1f%% Confidence Interval", x$level*100))
	if (is.na(x$L2)) {
#		out <- array(data=c(length(x$Y1),x$estimate, x$L1,x$U1), c(1,4))
#		colnames(out) <- c("n", "Lower", "d0_hat", "Upper")
#		cat(sprintf("n=%d samples, d0_hat=%.3f, CI: [%.3f,%.3f]\n",length(x$Y1),x$estimate, x$L1,x$U1))
		cat(sprintf("\nn   Lower   d0_hat   Upper\n%d   %.3f   %.3f   %.3f\n",length(x$Y1),x$L1,x$estimate,x$U1))
	}
	else {
		cat(sprintf("\nn1   n2   Lower   d0_hat   Upper\n%d   %d   %.3f   %.3f   %.3f\n",length(x$Y1),length(x$Y2),x$L2,x$estimate,x$U2))
#		cat(sprintf("n1=%d,n2=%d samples yields CI: [%.3f,%.3f], with d0_hat=%.3f\n",length(x$Y1), length(x$Y2),x$L2,x$U2, x$estimate))
	}
  	invisible(x)
}

