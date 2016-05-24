##' Convert cohen's d into r
##'
##' Convert cohen's d into r
##'	
##' @param d cohen's d
##' @param p the proportion of individuals in group 1
##' @seealso \code{\link{r.to.d}}
##' @author Dustin Fife
##' @export
##' @examples
##' d.to.r(1.1, p=.3)
d.to.r = function(d, p=.5){
	sqrt(p*(1-p))*d/sqrt(p*(1-p)*d^2 + 1)
}
##' Convert r into cohen's d
##'
##' Convert r into cohen's d
##'	
##' @param r the correlation coefficient
##' @param p the proportion of individuals in group 1
##' @seealso \code{\link{d.to.r}}
##' @author Dustin Fife
##' @export
##' @examples
##' r.to.d(.3, p=.3)
r.to.d = function(r, p=.5){
	r/sqrt((p*(1-p))*(1-r^2))
}