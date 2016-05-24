function.MinValuePPV <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction) 
	if (control$valuePPV < 0 || control$valuePPV > 1) {
		stop("You have entered an invalid minimum value for Positive Predictive Value. \n The minimum value for Positive Predictive Value must be between 0 and 1.", call. = FALSE)
	}
	if (control$valuePPV == 0) {
		  warning ("You have entered the minimum possible value for Positive Predictive \n Value. All the cutpoints fulfill the condition. Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	if (control$valuePPV == 1) {
 		warning ("You have entered the maximum possible value for Positive Predictive \n Value. Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	index.cutpointsPPV <- which(measures.acc$PPV[,1] >= control$valuePPV) 
	if (length(index.cutpointsPPV) == 0) {
		 warning("There is no cutoff that fulfills this condition. Please introduce another value, if desired.", call. = FALSE, immediate. = TRUE)
		cMinValuePPV <- NULL
	} 
	if (length(index.cutpointsPPV)!= 0) {
		cutpointsPPV <- measures.acc$cutoffs[index.cutpointsPPV]		
		if (length(index.cutpointsPPV) == 1) {
			cMinValuePPV <- cutpointsPPV
		}
		if (length(index.cutpointsPPV)> 1) {
			NPVnew <- obtain.optimal.measures(cutpointsPPV, measures.acc)$NPV				
			cutpointsNPVnew <- cutpointsPPV[which(round(NPVnew[,1],10) == round(max(NPVnew[,1], na.rm = TRUE),10))] 
					  
			if (length(cutpointsNPVnew)> 1) {
				PPVnew <- obtain.optimal.measures(cutpointsNPVnew, measures.acc)$PPV
				cMinValuePPV <- cutpointsNPVnew[which(round(PPVnew[,1],10) == round(max(PPVnew[,1], na.rm = TRUE),10))]									 
			}
			if (length(cutpointsNPVnew) == 1) {
					cMinValuePPV <- cutpointsNPVnew
			}
		}
	}
	optimal.cutoff <- obtain.optimal.measures(cMinValuePPV, measures.acc)  
	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
