# Copyright (C) 2011 Pierrick Bruneau, see README for full notice

displayGraph <- function(measure, dev, vect, xlab="K", ylab="measure", main=" ") {
	# displays the provided measures with variability in shaded style
	
	ylim <- c(min(measure -dev), max(measure+dev))
	plot(vect, measure,ylim=ylim,type="n", xlab=xlab, ylab=ylab, main=main)
	polygon( c(vect,rev(vect)),
	         c(measure-dev,rev(measure+dev)),
		 col = "gray90",
		 lty=2,border=NA)
	lines( vect , measure - dev , lty=2)
	lines( vect, measure+dev , lty=2)
	
	
	
	lines( vect, measure , lwd=2 , col="black")
	
}
