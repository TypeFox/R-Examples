##' Compute the unrestricted variance of T.
##'
##' Given reliability estimates in the population and in the sample, as well as the variance of X, \code{Ut} estimates
##' the unrestricted variance of T.
##'	
##' @param rxxi the reliability estimate of x in the incumbent (restricted) population
##' @param rxxa the reliability estimate of x in the applicant (unrestricted) population
##' @param ux The unrestricted variance of X
##' @seealso \code{\link{caseIV}}
##' @return a value indicating the unresticted variance of T
##' @author Dustin Fife
##' @export
##' @examples
##' Ut(.6, .8, 1.5)
Ut = function(rxxi, rxxa, ux){
	ut = sqrt(ux^2*(rxxi/rxxa))
	return(ut)
}