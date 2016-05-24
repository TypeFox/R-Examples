ordinDNA <- function(distobj, sppVector, ...){
	#Conduct Principal Coordinates Analysis
	pco <- cmdscale(distobj, length(sppVector) - 1, eig = TRUE)
	
	ordinDNA <- list(pco = pco, sppVector = sppVector)
	
	attr(ordinDNA, "class") <- "ordinDNA"
	
	#Plotting
	plot.ordinDNA(ordinDNA, ...)
	
	#Return object
	invisible(ordinDNA)
}
