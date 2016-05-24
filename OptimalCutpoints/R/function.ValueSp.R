function.ValueSp <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$valueSp < 0 || control$valueSp > 1) {
		stop("You have entered an invalid value for Specificity. \n The value for Specificity must be between 0 and 1.", call. = FALSE)		
	}
	if (control$valueSp == 0) {
		warning("You have entered the minimum possible value for Specificity. \n Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	if (control$valueSp == 1) {
		warning("You have entered the maximum possible value for Specificity. \n Please check this value.", call. = FALSE, immediate. = TRUE)
	}
	index.cutpoints <- which(round(measures.acc$Sp[,1],10) == round(control$valueSp,10))					  
	if (length(index.cutpoints)== 0) {		 
		warning("There is no cutpoint that yields the exact Specificity designated. The cutpoint having the closest value to the designated Specificity has therefore been selected.", call. = FALSE, immediate. = TRUE)			
		difference <- abs(control$valueSp-measures.acc$Sp[,1])
		index.cutpoints <- which(round(difference,10) == round(min(difference,na.rm=TRUE),10))							
	}
	if (length(index.cutpoints)!= 0) {
		if (length(index.cutpoints)== 1) {
			cvalueSp <- measures.acc$cutoffs[index.cutpoints]
		}
		if (length(index.cutpoints)!= 1) {
			cutpoints <- measures.acc$cutoffs[index.cutpoints]		  
			Senew <- obtain.optimal.measures(cutpoints, measures.acc)$Se		   
			cutpointsSenew <- cutpoints[which(round(Senew[,1],10) == round(max(Senew[,1],na.rm=TRUE),10))]				
			cvalueSp <- cutpointsSenew
		}   
	}

	optimal.cutoff <- obtain.optimal.measures(cvalueSp, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
