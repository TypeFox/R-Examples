summary.bayesQR <- function(object, burnin=0, credint=c(.025,.975), quantile=NULL, ...){

  # Create function for error message
  pandterm <- function(message) {
    stop(message, call. = FALSE)
  }

  # If no quantile is specified, summarize all estimated quantiles 
	if (is.null(quantile)){
    out <- lapply(object, FUN="summary.bayesQR.single", burnin=burnin, credint=credint)

	# Else, only summarize subset of estimated quantiles
	} else {

	  # Specified quantiles should be in object
	  if (!all(quantile %in% sapply(object, "[[", "quantile"))){
		  pandterm("One or more specified quantiles were not estimated")
		} else {

		  # Subset object
		  object <- object[(sapply(object, "[[", "quantile") %in% quantile)]

			# Summarize specified quantiles
      out <- lapply(object, FUN="summary.bayesQR.single", burnin=burnin, credint=credint)
		}
	}

  class(out) <- "bayesQR.summary"
  return(out)
}
