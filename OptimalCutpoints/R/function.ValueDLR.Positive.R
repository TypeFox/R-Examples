function.ValueDLR.Positive <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$valueDLR.Positive < 0) {
		stop("You have entered an invalid value for the Positive Diagnostic Likelihood \n Ratio. The Positive Diagnostic Likelihood Ratio must be positive.", call. = FALSE)
	}

	cValueDLR.Positive <- measures.acc$cutoffs[which(round(measures.acc$DLR.Positive[,1],10) == round(control$valueDLR.Positive,10))]  
 
	if (length(cValueDLR.Positive)== 0) {
		warning("There is no cutpoint that yields the exact Diagnostic Positive \n Likelihood Ratio designated. The cutpoint having the closest value to the \n designated Diagnostic Positive Likelihood Ratio has therefore been selected.", call. = FALSE, immediate. = TRUE)
		difference <- abs(control$valueDLR.Positive-measures.acc$DLR.Positive[,1])
		index.cutpoints <- which(round(difference,10) == round(min(difference, na.rm = TRUE),10)) 
 				  
		cValueDLR.Positive <- measures.acc$cutoffs[index.cutpoints]		 
 	}	

	optimal.cutoff <- obtain.optimal.measures(cValueDLR.Positive, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
