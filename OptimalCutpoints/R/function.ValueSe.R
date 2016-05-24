function.ValueSe <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$valueSe < 0 || control$valueSe > 1) {
		stop("You have entered an invalid value for Sensitivity. \n The value for Sensitivity must be between 0 and 1.", call. = FALSE)		
	}
	if (control$valueSe == 0) {
		warning("You have entered the minimum possible value for Sensitivity. \n Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	if (control$valueSe == 1) {
		warning("You have entered the maximum possible value for Sensitivity. \n Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	
	index.cutpoints <- which(round(measures.acc$Se[,1],10) == round(control$valueSe,10))  
	
	if (length(index.cutpoints) == 0) {		 
		warning("There is no cutpoint that yields the exact Sensitivity designated. \n The cutpoint having the closest value to the designated Sensitivity has therefore been selected.", call. = FALSE, immediate. = TRUE)
		difference <- abs(control$valueSe-measures.acc$Se[,1])
		index.cutpoints <- which(round(difference,10) == round(min(difference,na.rm=TRUE),10))							 
	}	 
	if (length(index.cutpoints)!= 0) {
		if (length(index.cutpoints)== 1) {
			cvalueSe <- measures.acc$cutoffs[index.cutpoints]
		}
		if (length(index.cutpoints)!= 1) {
			cutpoints <- measures.acc$cutoffs[index.cutpoints]
			Spnew <- obtain.optimal.measures(cutpoints, measures.acc)$Sp
			cutpointsSpnew <- cutpoints[which(round(Spnew[,1],10) == round(max(Spnew[,1],na.rm=TRUE),10))]
		 	cvalueSe <- cutpointsSpnew		
		}	
	}

	optimal.cutoff <- obtain.optimal.measures(cvalueSe, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
