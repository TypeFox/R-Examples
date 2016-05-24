#' @title Creates a scan.gwaa object
#'
#' @description
#' Creates a scan.gwaa object from a GenABEL object and p-values.
#'
#' @param data A gwaa.data object
#' @param P1df P-values computed from external analysis
#' @param SNP.eff Estimated additive SNP effects
#' 
#' @author Lars Ronnegard
#' 
Create_gwaa_scan <- function(data, P1df, SNP.eff) {
	chi2.1df <- qchisq(P1df, df=1, lower.tail=FALSE)
	lambda <- try( estlambda(chi2.1df), silent=TRUE )
	se_effB <- abs(SNP.eff)/qnorm(P1df/2, lower.tail=FALSE)
    if ( class(lambda) == "try-error" ) lambda <- list( estimate = NA, se = NA )
	results <- data.frame(effB = SNP.eff,  se_effB = se_effB, chi2.1df = NA, P1df = P1df, Pc1df = NA, stringsAsFactors = FALSE)
  rownames(results) = snpnames(data)
	out <- new("scan.gwaa", results = results, annotation = annotation(data), lambda = lambda, idnames = idnames(data), call = as.call( list("rGLS") ), family = "gaussian")
	return(out)
}

