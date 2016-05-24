function.ValueNPV <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$valueNPV < 0 || control$valueNPV > 1) {
		stop("You have entered an invalid value for Negative Predictive Value. \n The value for Negative Predictive Value must be between 0 and 1.", call. = FALSE)		
	}
	if (control$valueNPV == 0) {
		warning("You have entered the minimum possible value for Negative Predictive Value. \n Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	if (control$valueNPV == 1) {
		warning("You have entered the maximum possible value for Negative Predictive Value. \n Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	
	index.cutpoints <- which(round(measures.acc$NPV[,1],10) == round(control$valueNPV,10))					  
	if (length(index.cutpoints)== 0) {		 
		warning("There is no cutpoint that yields the exact Negative Predictive Value designated. The cutpoint having the closest value to the designated Negative Predictive Value has therefore been selected.", call. = FALSE, immediate. = TRUE)			
		difference <- abs(control$valueNPV-measures.acc$NPV[,1])
		index.cutpoints <- which(round(difference,10) == round(min(difference,na.rm=TRUE),10))							
	}
	if (length(index.cutpoints)!= 0) {
		if (length(index.cutpoints)== 1) {
			cvalueNPV <- measures.acc$cutoffs[index.cutpoints]
		}
		if (length(index.cutpoints)!= 1) {
			cutpoints <- measures.acc$cutoffs[index.cutpoints]		  
			PPVnew <- obtain.optimal.measures(cutpoints, measures.acc)$PPV		   
			cutpointsPPVnew <- cutpoints[which(round(PPVnew[,1],10) == round(max(PPVnew[,1],na.rm=TRUE),10))]				
			cvalueNPV <- cutpointsPPVnew
		}   
	}

	optimal.cutoff <- obtain.optimal.measures(cvalueNPV, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
