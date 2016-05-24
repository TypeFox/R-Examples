function.MinValueNPV <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$valueNPV < 0 || control$valueNPV > 1) {
		stop("You have entered an invalid minimum value for Negative Predictive Value. \n The minimum value for Negative Predictive Value must be between 0 and 1.", call. = FALSE)
	}
	if (control$valueNPV == 0) {
		warning ("You have entered the minimum possible value for Negative Predictive \n Value. All the cutpoints fulfill the condition. Please check this value.", call. = FALSE, immediate. = TRUE)	 		
	}
	if (control$valueNPV == 1) {
		warning ("You have entered the maximum possible value for Negative Predictive \n Value. Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	
	index.cutpointsNPV <- which(measures.acc$NPV[,1] >= control$valueNPV)		 
	if (length(index.cutpointsNPV) == 0) {
		warning("There is no cutoff that fulfills this condition. Please introduce another value, if desired.", call. = FALSE, immediate. = TRUE)
		cMinValueNPV <- NULL
	} 
	if (length(index.cutpointsNPV)!= 0) {
		cutpointsNPV <- measures.acc$cutoffs[index.cutpointsNPV]
	
		if (length(index.cutpointsNPV) == 1) {
			cMinValueNPV <- cutpointsNPV
		}
		if (length(index.cutpointsNPV)> 1) {
			PPVnew <- obtain.optimal.measures(cutpointsNPV, measures.acc)$PPV 		 		
			cutpointsPPVnew <- cutpointsNPV[which(round(PPVnew[,1],10) == round(max(PPVnew[,1], na.rm = TRUE),10))] 
					  
			if (length(cutpointsPPVnew)> 1) {
					NPVnew <- obtain.optimal.measures(cutpointsPPVnew, measures.acc)$NPV
					cMinValueNPV <- cutpointsPPVnew[which(round(NPVnew[,1],10) == round(max(NPVnew[,1], na.rm = TRUE),10))]				 		
			}
			if (length(cutpointsPPVnew) == 1) {
					cMinValueNPV <- cutpointsPPVnew
			}
		}
  	}
	optimal.cutoff <- obtain.optimal.measures(cMinValueNPV, measures.acc)
	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
