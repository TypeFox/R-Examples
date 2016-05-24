#' Estimate correlations with the EM
#' 
#' This function estimates unattenuated correlations using the EM algorithm.
#' 
#' @param data.matrix a matrix where the rows consist of observations and the
#' columns consist of variables
#' @return returns a corrected correlation between the variables
#' @author Dustin Fife
#' @seealso See Also as \code{\link{rel.correction}} to correct the estimate
#' for unreliability, \code{\link{caseIV}}, \code{\link{caseIIIR}},
#' \code{\link{caseIII}}
#' @export
#' @import norm
#' @examples
#' 
#' data(selection.example.data)
#' corrected.cor.matrix = em(data.matrix(selection.example.data))
#' corrected.cor.matrix
em <-
function(data.matrix){
	ss = prelim.norm(data.matrix)
	thetahat = em.norm(ss, showits=FALSE)
	cor.corrected = data.frame(getparam.norm(ss, thetahat, corr=TRUE)$r)
	if (!is.null(names(data.matrix))){
		names(cor.corrected) = names(data.matrix)
		row.names(cor.corrected) = names(data.matrix)
	}

	cat("Returning corrected correlation matrix:\n")
	return(cor.corrected)
}
