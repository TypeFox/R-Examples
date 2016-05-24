function.MaxDOR <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)		
	TP <- measures.acc$Se[,1]*measures.acc$n$d
	TP[TP == 0] <- 0.5
	
	TN <- measures.acc$Sp[,1]*measures.acc$n$h
	TN[TN == 0] <- 0.5
	
	FN <- (1-measures.acc$Se[,1])*measures.acc$n$d
	FN[FN == 0] <- 0.5
	
	FP <- (1-measures.acc$Sp[,1])*measures.acc$n$h
	FP[FP == 0] <- 0.5
	
	DOR <- (TP*TN)/(FN*FP)

	cMaxDOR <- measures.acc$cutoffs[which(round(DOR,10) == round(max(DOR, na.rm = TRUE),10))]
	optimal.DOR <- max(DOR, na.rm = TRUE)

	optimal.cutoff <- obtain.optimal.measures(cMaxDOR, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = DOR, optimal.criterion = optimal.DOR)
	res
}
