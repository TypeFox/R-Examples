function.MinValueSp <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$valueSp < 0 || control$valueSp > 1) {
		stop("You have entered an invalid minimum value for Specificity. \n The minimum value for Specificity must be between 0 and 1.", call. = FALSE)		
	}
	if (control$valueSp == 0) {
		warning("You have entered the minimum possible value for Specificity. \n All the cutpoints fulfill the condition. Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	if (control$valueSp == 1) {
		warning("You have entered the maximum possible value for Specificity. \n Please check this value.", call. = FALSE, immediate. = TRUE)		
	}
	index.cutpoints <- which(measures.acc$Sp[,1] >= control$valueSp)
	if (length(index.cutpoints) == 0)
	{
	warning("There is no cutoff that fulfills this condition. Please, enter a new value, if desired.", call. = FALSE, immediate. = TRUE)
	cMinValueSp <- NULL		
	}
	if (length(index.cutpoints)!= 0) {
		cutpoints <- measures.acc$cutoffs[index.cutpoints]		  
		if (length(index.cutpoints) == 1) {
			cMinValueSp <- cutpoints
		}
		if (length(index.cutpoints)!= 1) {
			Senew <- obtain.optimal.measures(cutpoints, measures.acc)$Se
			cutpointsSenew <- cutpoints[which(round(Senew[,1],10) == round(max(Senew[,1],na.rm=TRUE),10))] 
															   
			if (length(cutpointsSenew)> 1) {
				Spnew <- obtain.optimal.measures(cutpointsSenew, measures.acc)$Sp
				cMinValueSp <- cutpointsSenew[which(round(Spnew[,1],10) == round(max(Spnew[,1],na.rm=TRUE),10))]	
			}	 
			if (length(cutpointsSenew)== 1) {
				cMinValueSp <- cutpointsSenew
			}
		}
	}
	optimal.cutoff <- obtain.optimal.measures(cMinValueSp, measures.acc)
	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
