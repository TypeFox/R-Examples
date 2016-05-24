#' Correct x/y correlations using Case IV
#' 
#' \code{caseIV} uses Hunter and Schmidt's (2004) correction procedure to
#' estimate the correlation between t and p, assuming direct selection on a
#' latent s.
#' 
#' @param rtpi the correlation between t and p in the incumbent pool.
#' @param rxxi the reliability of x in the incumbent pool.
#' @param rxxa the reliability of x in the applicant pool.
#' @param ux the ratio of selected to unselected variance in x.
#' @return a scalar that is the estimate of the correlation between t and p.
#' @author Dustin Fife
#' @export
#' @seealso See Also \code{\link{rel.correction}} to obtain estimates of rtpi
#' from rxyi, \code{\link{caseIII}}, \code{\link{em}}, \code{\link{caseIIIR}}
#' @references Hunter, J.E., and Schmidt, F.L. (2004). Methods of
#' meta-analysis: Correcting error and bias in research findings. Thousand
#' Oaks, CA: Sage.
#' @examples
#' 
#' 	# correct an xy correlation for unreliability
#' 	rtpi = rel.correction(rxyi=.3, rxxi=.8, ryyi=.6)
#' 
#' 	# estimate Case IV
#' 	caseIV(rtpi=rtpi, rxxi=.6, rxxa=.8, ux=.8)
caseIV <-
function(rtpi,rxxi, rxxa, ux){
	ut = 1/Ut(rxxi, rxxa, ux)
	c.IV = (ut*rtpi)/sqrt(1 + (ut^2-1)*rtpi^2)
	return(c.IV)
}