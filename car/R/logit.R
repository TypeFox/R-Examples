# logit transformation of proportion or percent (J. Fox)

# last modified 2012-06-24 by J. Fox

logit <- function(p, percents=range.p[2] > 1, adjust){
	range.p <- range(p, na.rm=TRUE)
	if (percents){
		if (range.p[1] < 0 || range.p[1] > 100) stop("p must be in the range 0 to 100")
		p <- p/100
		range.p <- range.p/100
	}
	else if (range.p[1] < 0 || range.p[1] > 1) stop("p must be in the range 0 to 1")
	a <-if (missing(adjust)) {
				if (isTRUE(all.equal(range.p[1], 0)) || isTRUE(all.equal(range.p[2], 1))) .025 else 0
			}
			else adjust
	if (missing(adjust) && a != 0) warning(paste("proportions remapped to (", a, ", ", 1-a, ")", sep=""))
	a <- 1 - 2*a
	log((0.50 + a*(p - 0.50))/(1 - (0.50 + a*(p - 0.50))))
}
    
