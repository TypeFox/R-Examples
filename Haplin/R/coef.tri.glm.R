coef.tri.glm <- function(object, cov.type = "Fisher.EM", ...){
##
## EXTRACTS COEFFICIENTS AND THE VARIANCE-COVARIANCE MATRIX FROM A tri.glm OBJECT
##
#
## EXTRACT glm PART
.res <- object$result
#
.summary <- summary.glm(.res, correlation = F)
#
## COEFFICIENTS
.coef <- .summary$coefficients[, 1] # WORKS FOR BOTH SPLUS AND R
#
## HIERARCHICAL SELECTION OF COVARIANCE MATRIX
if(cov.type == "resamp"){
	## RESAMPLED (JACKKNIFE) MATRIX
	.cov <- attr(object, "cov.resamp")
	if(is.null(.cov)){
		warning("Resampled (jackknife) covariance matrix not available!\n Using standard Fisher EM-corrected instead")
		cov.type <- "Fisher.EM"
	}
}
if(cov.type == "Fisher.EM"){
	## FISHER COVARIANCE MATRIX, CORRECTED FOR EM-IMPUTATIONS
	.cov <- attr(object, "cov.correct")
	if(is.null(.cov)){
		warning("Fisher EM-corrected covariance not available!\n Using uncorrected instead")
		cov.type <- "Fisher"
	}
}
if(cov.type == "Fisher"){
	## STANDARD GLM FISHER COVARIANCE MATRIX, NOT CORRECTED FOR EM-IMPUTATIONS
	.cov <- .summary$cov.unscaled	##
}
#
##
attr(.cov, "cov.type") <- cov.type
.ut <- list(coef = .coef, cov = .cov)
#
return(.ut)
}
