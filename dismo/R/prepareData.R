# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : June 2009
# Version 1
# Licence GPL v3

prepareData <- function(x, p, b, factors, xy=FALSE) {
	if (xy) {
		coords <- data.frame(rbind(as.matrix(p), as.matrix(b)))
		colnames(coords) <- c('x', 'y')
	} 

	p <- extract(x, p)
	b <- extract(x, b)
	pb <- data.frame(  cbind(pb=c(rep(1, nrow(p)), rep(0, nrow(b))), rbind(p, b)) )
	if (!missing(factors)) {
		for (f in factors) {
			pb[,f] = factor(pb[,f])
		}
	}
	
	if (xy) {
		pb <- cbind(coords, pb)
	}
	return(pb)
}

