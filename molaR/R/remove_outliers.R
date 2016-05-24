#' Mask outliers on some faces
#'
#' This function will block out the top 0.1 percent of the faces
#' @param Energy_values energy density values on faces
#' @param X percentile above which to remove
#'
#' 
#' remove_outliers()

remove_outliers <- function(Energy_values, X) {
	DNEs <- Energy_values$DNE_Values
	FAs <- Energy_values$Face_Areas
	Q <- quantile(DNEs, probs=c(X))
	
	OutlierList <- which(DNEs > Q)
	
	DNEs[OutlierList] <- 0
	
	out <- data.frame(DNE_Values=DNEs, Face_Areas=FAs)
	
	return(out)
}	
