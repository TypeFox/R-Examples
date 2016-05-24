#' Plot generic for SmithWilsonYieldCurve objects
#' 
#' @param x An object of class SmithWilsonYieldCurve or a vector of terms to evaluate the curve at
#' @param y Optionally an object of class SmithWilsonYieldCurve
#' @param aspect either "cts" for continously compounded spot rates, or "zero" for ZCB prices
#' @param ... other arguments to pass to the default plot function
#' 
#' @method plot SmithWilsonYieldCurve
#' @export
#' 
plot.SmithWilsonYieldCurve <- function(x, y, ..., aspect=c("cts", "zero")){
	
	if (missing(y)) 
		terms <- 1:50
	 else 
		terms <- y
		
	zeros <- x$P(terms)
	
	switch( match.arg(aspect),
			cts = plot( terms, -log(zeros) / terms, xlab="Term", ylab="Rate", ... ),
			zero = plot( terms, zeros, xlab="Term", ylab="ZCB price", ... ),
			stop("Unknown aspect of yield curve to plot")
	)
	
}