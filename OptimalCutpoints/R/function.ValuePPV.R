function.ValuePPV <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$valuePPV < 0 || control$valuePPV > 1) {
		stop("You have entered an invalid value for Positive Predictive Value. \n The value for Positive Predictive Value must be between 0 and 1.", call. = FALSE)		
	}
	if (control$valuePPV == 0) {
		warning("You have entered the minimum possible value for Positive Predictive Value. \n Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	if (control$valuePPV == 1) {
		warning("You have entered the maximum possible value for Positive Predictive Value. \n Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	
	index.cutpoints <- which(round(measures.acc$PPV[,1],10) == round(control$valuePPV,10))  
	
	if (length(index.cutpoints) == 0) {		 
		warning("There is no cutpoint that yields the exact Positive Predictive Value designated. \n The cutpoint having the closest value to the designated Positive Predictive Value has therefore been selected.", call. = FALSE, immediate. = TRUE)
		difference <- abs(control$valuePPV-measures.acc$PPV[,1])
		index.cutpoints <- which(round(difference,10) == round(min(difference,na.rm=TRUE),10))							 
	}	 
	if (length(index.cutpoints)!= 0) {
		if (length(index.cutpoints)== 1) {
			cvaluePPV <- measures.acc$cutoffs[index.cutpoints]
		}
		if (length(index.cutpoints)!= 1) {
			cutpoints <- measures.acc$cutoffs[index.cutpoints]
			NPVnew <- obtain.optimal.measures(cutpoints, measures.acc)$NPV
			cutpointsNPVnew <- cutpoints[which(round(NPVnew[,1],10) == round(max(NPVnew[,1],na.rm=TRUE),10))]
		 	cvaluePPV <- cutpointsNPVnew		
		}	
	}

	optimal.cutoff <- obtain.optimal.measures(cvaluePPV, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
