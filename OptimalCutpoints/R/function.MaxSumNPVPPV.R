function.MaxSumNPVPPV <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)	
	sum <- measures.acc$PPV[,1] + measures.acc$NPV[,1]	
	cmaxSumNPVPPV <- measures.acc$cutoffs[which(round(sum,10) == round(max(sum,na.rm=TRUE),10))]
	optimal.sum <- max(sum,na.rm=TRUE)

	optimal.cutoff <- obtain.optimal.measures(cmaxSumNPVPPV, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = sum, optimal.criterion = optimal.sum)
	res
}
