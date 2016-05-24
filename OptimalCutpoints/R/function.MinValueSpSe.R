function.MinValueSpSe <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$valueSp < 0 || control$valueSp > 1) {
		stop("The minimum value for Specificity must be between 0 and 1.", call. = FALSE)
	}
	if (control$valueSe < 0 || control$valueSe > 1) {
		stop("The minimum value for Sensitivity must be between 0 and 1.", call. = FALSE)
	}
	if (control$valueSp == 0 & control$valueSe == 0) {
		warning ("You have entered the minimum possible values for Specificity and \n Sensitivity. All the cutpoints fulfill the condition. Please check these values.", call. = FALSE, immediate. = TRUE)
	}
	if (control$valueSp == 1 & control$valueSe == 1) {
		warning ("You have entered the maximum possible values for Specificity and \n Sensitivity. Please check these values.", call. = FALSE, immediate. = TRUE)
	}
	if (is.logical(control$maxSp) == FALSE) {
		stop("'maxSp' must be a logical-type argument.", call. = FALSE)		 	
	}
	index.cutpoints <- which(measures.acc$Sp[,1] >= control$valueSp & measures.acc$Se[,1] >= control$valueSe)	
	if (length(index.cutpoints)== 0) {
	 	warning("There is no cutoff that fulfills these conditions. \n Please enter other minimum values, if desired.", call. = FALSE, immediate. = TRUE)
		cMinValueSpSe <- NULL
	}
	if (length(index.cutpoints)!= 0) {
		if (length(index.cutpoints)== 1) {
			cMinValueSpSe <- measures.acc$cutoffs[index.cutpoints]
		}
		# If there is more than one cutpoint fulfilling these conditions, those which yield 
		# maximum Sensitivity or maximum Specificity are chosen:
		if (length(index.cutpoints)> 1) {
			cutpoints <- measures.acc$cutoffs[index.cutpoints]	
			### If you seek to maximize Specificity:			 	
			if(control$maxSp == TRUE) {
				Spnew <- obtain.optimal.measures(cutpoints, measures.acc)$Sp
				cutpointsSpnew <- cutpoints[which(round(Spnew[,1],10) == round(max(Spnew[,1],na.rm=TRUE),10))] 
							  
				if (length(cutpointsSpnew)> 1) {
					Senew <- obtain.optimal.measures(cutpointsSpnew, measures.acc)$Se
					cMinValueSpSe <- cutpointsSpnew[which(round(Senew[,1],10) == round(max(Senew[,1],na.rm=TRUE),10))]					 		
				}
				if (length(cutpointsSpnew)== 1) {
					cMinValueSpSe <- cutpointsSpnew
				}
			}
			### If you seek to maximize Sensitivity:
			if(control$maxSp == FALSE) {
				Senew <- obtain.optimal.measures(cutpoints, measures.acc)$Se
				cutpointsSenew <- cutpoints[which(round(Senew[,1],10) == round(max(Senew[,1],na.rm=TRUE),10))] 
								
				if (length(cutpointsSenew)> 1) {
					Spnew <- obtain.optimal.measures(cutpointsSenew, measures.acc)$Sp
					cMinValueSpSe <- cutpointsSenew[which(round(Spnew[,1],10) == round(max(Spnew[,1],na.rm=TRUE),10))]					 		
				}
				if (length(cutpointsSenew)== 1) {
					cMinValueSpSe <- cutpointsSenew
				}
			}
		}
	}
	optimal.cutoff <- obtain.optimal.measures(cMinValueSpSe, measures.acc)
	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
