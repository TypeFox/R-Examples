bayesQR <- function(formula=NULL, data=NULL, quantile=0.5, alasso=FALSE, ndraw=NULL, keep=1, prior=NULL){

  # Create function for error message
  pandterm <- function(message) {
    stop(message, call. = FALSE)
  }

  # Check if a single or a sequence of quantile regressions is required
	nqr <- length(quantile)

  # Define an empty object (list)
	out <- NULL

	# If only one quantile is required, then call bayesQR.single
	if (nqr==1){
	  out[[1]] <- bayesQR.single(formula=formula, data=data, quantile=quantile, alasso=alasso, ndraw=ndraw, keep=keep, prior=prior)

	# Else, estimate a sequence of bayesQR.single
	} else {

		# Sort required quantiles
	  quantile <- sort(quantile)

		# Estimate sequence of bayesQR
		for (i in 1:nqr){

      # Print information to console
      cat("************************************************","\n")
      cat("* Start estimating quantile ", i," of ", nqr, "in total *", "\n")
      cat("************************************************","\n")
  
      # Set correct quantile and estimate model
	    out[[i]] <- bayesQR.single(formula=formula, data=data, quantile=quantile[i], alasso=alasso, ndraw=ndraw, keep=keep, prior=prior)
		}
	}
  class(out) <- "bayesQR"
	return(out)
}
