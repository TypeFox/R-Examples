#' Means between change-points
#'
#' The function finds the average of the input vector \code{x} between change-points given in \code{cpt}.    
#'
#' @param x a vector
#' @param cpt a vector of integers with localisations of change-points
#' @param ... further arguments passed to \code{mean} method
#' @return a vector of the same length as \code{x}, piecewise constant and equal to the mean between change-points given in \code{cpt} 
#' @export means.between.cpt
#' @examples
#' x <- rnorm(100)+c(rep(-1,50),rep(1,50))
#' cpt <- 50
#' means.between.cpt(x,cpt)
#' w <- wbs(x)
#' cpt <- changepoints(w)
#' means.between.cpt(x,cpt=cpt$cpt.ic$sbic)

means.between.cpt <-
		function(x, cpt=NULL,...) {
	
	x <- as.numeric(x)
	if(NA%in%x) stop("x vector cannot contain NA's")
	
	
	n <- length(x)
	
	if(!is.null(cpt)) cpt <- as.integer(cpt)
	cpt <- cpt[!is.na(cpt)]
	
	
	
	len.cpt <- length(cpt)
	
	
	if (len.cpt) {
		if(min(cpt)<1 || max(cpt)>=n) stop("change-points should be between 1 and n-1")
		cpt <- sort(cpt)
	}
	
	s <- e <- rep(0, len.cpt+1)
	s[1] <- 1
	e[len.cpt+1] <- n
	if (len.cpt) {
		s[2:(len.cpt+1)] <- cpt+1
		e[1:len.cpt] <- cpt
	}
	means <- rep(0, len.cpt+1)
	for (i in 1:(len.cpt+1)) means[i] <- mean(x[s[i]:e[i]],...)
	rep(means, e-s+1)
}