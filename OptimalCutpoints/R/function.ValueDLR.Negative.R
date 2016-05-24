function.ValueDLR.Negative <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$valueDLR.Negative < 0) {
		stop("You have entered an invalid value for the Negative Diagnostic Likelihood \n Ratio. The Negative Diagnostic Likelihood Ratio must be positive.", call. = FALSE)			
	}
	
	cValueDLR.Negative <- measures.acc$cutoffs[which(round(measures.acc$DLR.Negative[,1],10) == round(control$valueDLR.Negative,10))]	

	if (length(cValueDLR.Negative)== 0) {
		warning("There is no cutpoint that yields the exact Diagnostic Negative \n Likelihood Ratio designated. The cutpoint having the closest value to the \n designated Diagnostic Negative Likelihood Ratio has therefore been selected.", call. = FALSE, immediate. = TRUE)		
		difference <- abs(control$valueDLR.Negative-measures.acc$DLR.Negative[,1])
		index.cutpoints <- which(round(difference,10) == round(min(difference, na.rm = TRUE),10))
		cValueDLR.Negative <- measures.acc$cutoffs[index.cutpoints]
	}
	
	optimal.cutoff <- obtain.optimal.measures(cValueDLR.Negative, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
