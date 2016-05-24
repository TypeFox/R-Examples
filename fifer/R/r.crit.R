##' Compute the critical r value, or return the p value of an r, assuming a given number of degrees of freedom.
##'
##' @title Compute critical r or p.
##' @param n The sample size.
##' @param p the probability. Defaults to .025.
##' @param r The observed r value
##' @export
##' @return the critical r value or the observed p value for a given r
##' @author Dustin Fife
##' @examples
##' r.crit(n=100, p=.025)
##' r.crit(n=20, r=.6, p=.05)
r.crit = function(n, r=NULL, p=.025){
	
	if (is.null(r)){
		tcrit = qt(p, n)
		rcrit = sqrt(tcrit^2 / (n-2+tcrit^2))
		return(rcrit)
	} else {
		tobs = r*sqrt((n-2)/(1-r^2))
		pobt = 1-pt(tobs, df=n-2)
		tcrit = qt(p, n)
		rcrit = sqrt(tcrit^2 / (n-2+tcrit^2))		
		list(rcrit=rcrit, pobt=pobt)
	}
}
