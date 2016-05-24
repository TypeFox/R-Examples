##' Scale a variable to have a particular mean/sd (or min/max)
##'
##' @title Scale a variable
##' @param x the variable to be scaled
##' @param mean the mean you wish the distribution to have
##' @param sd the sd you wish the distribution to have
##' @param min the min you wish the distribution to have
##' @param max the max you wish the distribution to have
##' @return the scaled variable
##' @export
##' @author Dustin Fife
scaleIt = function(x, mean=NULL, sd=NULL, min=NULL, max=NULL){
	if (!is.null(min)){
		old.max = max(x, na.rm=T); old.min=min(x, na.rm=T)
		newx = (max-min)/(old.max-old.min)*(x-old.min)+min
		newx
	} else {
		newx = mean + sd*scale(x)
	}
}
