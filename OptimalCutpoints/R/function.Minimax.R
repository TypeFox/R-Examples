function.Minimax <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (is.logical(control$maxSp) == FALSE) {
		stop("'maxSp' must be a logical-type argument.", call. = FALSE)
	}
	FN <-(1-measures.acc$Se[,1])*length(data[data[,status] != tag.healthy, marker])
	FP <-(1-measures.acc$Sp[,1])*length(data[data[,status] == tag.healthy, marker])

	M <- vector()
	for(i in 1:length(measures.acc$cutoffs)) {
		if (FN[i] > FP[i]) {
			M[i] <- FN[i]
		} else  {
			M[i] <- FP[i]
		}
	}
	cMinimax <- measures.acc$cutoffs[which(round(M,10) == round(min(M,na.rm=TRUE),10))] 
	# If there is more than one cutpoint fulfilling these conditions, 
 	# those which yield maximum Sensitivity or maximum Specificity are chosen:
	if (length(cMinimax)> 1) {		
		### If you seek to maximize Specificity:
		if(control$maxSp == TRUE) {
			Spnew <- obtain.optimal.measures(cMinimax, measures.acc)$Sp			
			cutpointsSpnew <- cMinimax[which(round(Spnew[,1],10) == round(max(Spnew[,1],na.rm=TRUE),10))]
								  
			if (length(cutpointsSpnew)> 1) {
				Senew <- obtain.optimal.measures(cutpointsSpnew, measures.acc)$Se			 
				cMinimax <- cutpointsSpnew[which(round(Senew[,1],10) == round(max(Senew[,1],na.rm=TRUE),10))]			 		
			}
			if (length(cutpointsSpnew)== 1) {
				cMinimax <- cutpointsSpnew
			}
		}
		### If you seek to maximize Sensitivity:
		if(control$maxSp == FALSE) {
		 	Senew <- obtain.optimal.measures(cMinimax, measures.acc)$Se		 
		 	cutpointsSenew <- cMinimax[which(round(Senew[,1],10) == round(max(Senew[,1],na.rm=TRUE),10))]
					
		 	if (length(cutpointsSenew)> 1) {
				Spnew <- obtain.optimal.measures(cutpointsSenew, measures.acc)$Sp			 
				cMinimax <- cutpointsSenew[which(round(Spnew[,1],10) == round(max(Spnew[,1],na.rm=TRUE),10))]								 
		 	}
		 	if (length(cutpointsSenew)== 1) {
				cMinimax <- cutpointsSenew
		 	}
		}	 
	} 
	optimal.M <- min(M,na.rm=TRUE)   
		
	optimal.cutoff <- obtain.optimal.measures(cMinimax, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = M, optimal.criterion = optimal.M)
	res
}
