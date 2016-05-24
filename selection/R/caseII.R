#' Corrects correlations using Case II
#' 
#' Using Thorndike's Case II correction, \code{caseIII} corrects the xy
#' correlation for direct restriction on x
#' 
#' The Case II correction is defined as follows insert later The result is an
#' unbiased estimate of the unattenuated correlation between X and Y
#' 
#' @param rxy the restricted correlation between x (the indirectly selected
#' variable) and y (the outcome variable).
#' @param sx the unrestricted standard deviation of x
#' @param sxs the restricted standard deviation of x
#' @return a scalar that is the unbiased estimate of the correlation between X
#' and Y.
#' @author Dustin Fife
#' @export
#' @seealso \code{\link{caseIV}}, \code{\link{caseIIIR}}, \code{\link{em}},
#' \code{\link{rel.correction}}
#' @references Thorndike, R. L. (1949). Personnel selection: Test and
#' measurement techniques. Oxford, England: Wiley.
#' 
#' Pearson, K. (1903). Mathematical contributions to the theory of evolution.
#' XI. On the influence of natural selection on the variability and correlation
#' of organs. Philosophical Transactions of the Royal Society of London. Series
#' A, Containing Papers of a Mathematical or Physical Character, 200, 1-66.
#' @examples
#' 
#' 	# correct assuming direct selection on X
#' corrected = caseII(rxy=.5, sx=1.5, sxs=.9)	
#' corrected
#' 
caseII <-
function(rxy, sx, sxs){
	num = rxy*(sx/sxs)
	denom = sqrt(1-rxy^2 + rxy^2*(sx^2/sxs^2))
	return(num/denom)
}
