function.MinValueNPVPPV <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (is.logical(control$maxNPV) == FALSE) {
		stop("'maxNPV' must be a logical-type argument.", call. = FALSE)
	}
	if (control$valueNPV < 0 || control$valueNPV > 1) {
		stop("You have entered an invalid minimum value for Negative Predictive Value. \n The minimum value for Negative Predictive Value must be between 0 and 1.", call. = FALSE)
	} 
	if (control$valuePPV < 0 || control$valuePPV > 1) {
		stop("You have entered an invalid minimum value for Positive Predictive Value. \n The minimum value for Positive Predictive Value must be between 0 and 1.", call. = FALSE)
	}
	if (control$valueNPV == 0 & control$valuePPV == 0) {
	warning ("You have entered the minimum possible values for Predictive Values. \n All the cutpoints fulfill the condition. Please check these values.", call. = FALSE, immediate. = TRUE)			
	}
	if (control$valueNPV == 1 & control$valuePPV == 1) {
	warning ("You have entered the maximum possible values for Predictive Values. \n Please check these values.", call. = FALSE, immediate. = TRUE)
	}
	index.cutpoints <- which(measures.acc$PPV[,1] >= control$valuePPV & measures.acc$NPV[,1] >= control$valueNPV)	
	if (length(index.cutpoints) == 0) {
		warning("There is no cutoff that fulfills these conditions. Please introduce other values, if desired.", call. = FALSE, immediate. = TRUE)
		cMinValueNPVPPV <- NULL		
	}
	if (length(index.cutpoints)!= 0) {
		if (length(index.cutpoints) == 1) {
			cMinValueNPVPPV <- measures.acc$cutoffs[index.cutpoints]
		}
		if (length(index.cutpoints)> 1) {
			cutpoints <- measures.acc$cutoffs[index.cutpoints]				  
			### If you seek to maximize Negative Predictive Value:
			if (control$maxNPV == TRUE) {				 
				NPVnew <- obtain.optimal.measures(cutpoints, measures.acc)$NPV				
				cutpointsNPVnew <- cutpoints[which(round(NPVnew[,1],10) == round(max(NPVnew[,1], na.rm = TRUE),10))] 			   
				if (length(cutpointsNPVnew)> 1) {
					PPVnew2 <- obtain.optimal.measures(cutpointsNPVnew, measures.acc)$PPV
					cMinValueNPVPPV <- cutpointsNPVnew[which(round(PPVnew2[,1],10) == round(max(PPVnew2[,1], na.rm = TRUE),10))]								
				}
				if (length(cutpointsNPVnew)== 1) {
					cMinValueNPVPPV <- cutpointsNPVnew
				}
			}
			### If you seek to maximize Positive Predictive Value:				
			if (control$maxNPV == FALSE) {
				PPVnew <- obtain.optimal.measures(cutpoints, measures.acc)$PPV				
				cutpointsPPVnew <- cutpoints[which(round(PPVnew[,1],10) == round(max(PPVnew[,1], na.rm = TRUE),10))]
				if (length(cutpointsPPVnew)> 1) {
					NPVnew2 <- obtain.optimal.measures(cutpointsPPVnew, measures.acc)$NPV
					cMinValueNPVPPV <- cutpointsPPVnew[which(round(NPVnew2[,1],10) == round(max(NPVnew2[,1], na.rm = TRUE),10))]							
				}
				if (length(cutpointsPPVnew)== 1) {
					cMinValueNPVPPV <- cutpointsPPVnew
				}
			}
		}
	}
	optimal.cutoff <- obtain.optimal.measures(cMinValueNPVPPV, measures.acc)
	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
