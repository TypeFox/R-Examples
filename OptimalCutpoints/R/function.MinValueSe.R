function.MinValueSe <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)  
	if (control$valueSe < 0 || control$valueSe > 1) {
		stop("You have entered an invalid minimum value for Sensitivity. \n The minimum value for Sensitivity must be between 0 and 1.", call. = FALSE)
	}
	if (control$valueSe == 0) {
		warning("You have entered the minimum possible value for Sensitivity. \n All the cutpoints fulfill the condition. Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	if (control$valueSe == 1) {
		warning("You have entered the maximum possible value for Sensitivity. \n Please check this value.", call. = FALSE, immediate. = TRUE) 
	}		
	index.cutpoints <- which(measures.acc$Se[,1] >= control$valueSe)  
	if (length(index.cutpoints)== 0) {
		warning("There is no cutoff that fulfills this condition. Please, enter a new value, if desired.", call. = FALSE, immediate. = TRUE)
		cMinValueSe <- NULL
	}
	if (length(index.cutpoints)!= 0) {
		cutpoints <- measures.acc$cutoffs[index.cutpoints]		
		if (length(index.cutpoints)== 1) {
			cMinValueSe <- cutpoints
		}
		if (length(index.cutpoints)!= 1) {					 
			Spnew <- obtain.optimal.measures(cutpoints, measures.acc)$Sp
			cutpointsSpnew <- cutpoints[which(round(Spnew[,1],10) == round(max(Spnew[,1],na.rm=TRUE),10))]			
		 
			if (length(cutpointsSpnew)> 1) {
				Senew <- obtain.optimal.measures(cutpointsSpnew, measures.acc)$Se
				cMinValueSe <- cutpointsSpnew[which(round(Senew[,1],10) == round(max(Senew[,1],na.rm=TRUE),10))]				 		
			}

			if (length(cutpointsSpnew)== 1)  {
					cMinValueSe <- cutpointsSpnew
			}
		}
	}
	optimal.cutoff <- obtain.optimal.measures(cMinValueSe, measures.acc)
	
 	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
 	res
}
