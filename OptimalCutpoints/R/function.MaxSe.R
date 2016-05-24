function.MaxSe <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)   
	
	cutpointsSe <- measures.acc$cutoffs[which(round(measures.acc$Se[,1],10) == round(max(measures.acc$Se[,1],na.rm=TRUE),10))]
		 
	if (length(cutpointsSe)> 1) {
		Spnew <- obtain.optimal.measures(cutpointsSe, measures.acc)$Sp
		cMaxSe <- cutpointsSe[which(round(Spnew[,1],10) == round(max(Spnew[,1],na.rm=TRUE),10))]		
	}
	if (length(cutpointsSe)== 1) {
		cMaxSe <- cutpointsSe
	}

	optimal.cutoff <- obtain.optimal.measures(cMaxSe, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)
	res
}
