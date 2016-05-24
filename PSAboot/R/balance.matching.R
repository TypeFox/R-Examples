#' Returns balance for each covariate from propensity score matching.
#' 
#' @param index.treated a vector with the index of treated rows in \code{covs}.
#' @param index.control a vector with the index of control rows in \code{covs}.
#' @param covs data frame or matrix of covariates. Factors should already be recoded.
#'        See \code{\link{cv.trans.psa}}
#' @return a named vector with one element per covariate.
#' @export
balance.matching <- function(index.treated, index.control, covs) {
	if(length(index.treated) != length(index.control) & 
	   length(index.control) != nrow(covs)) {
		stop('The length of index.treated and index.control must be the same and equal to nrow(covs)!')
	}
	bal <- c()
	for(covar in names(covs)) {
		cov <- data.frame(Treated=covs[index.treated,covar],
						  Control=covs[index.control,covar])
		ttest <- t.test(cov$Treated, cov$Control, paired=TRUE)
		bal[covar] <- ttest$estimate / sd(c(cov[,1],cov[,2]))	
	}
	return(bal)
}
