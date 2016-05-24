function.ObservedPrev <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)	
	if (measures.acc$cutoffs < 0 || measures.acc$cutoffs > 1) {
		warning("Diagnostic marker values are not between 0 and 1 for this \n criterion. A data transformation has been performed.", call. = FALSE, immediate. = TRUE)
		tcutoffs <- (measures.acc$cutoffs- min(measures.acc$cutoffs))/(max(measures.acc$cutoffs)-min(measures.acc$cutoffs))
		difference <- abs(tcutoffs-calculate.sample.prev(data, status, tag.healthy))
	} else {
		difference <- abs(measures.acc$cutoffs-calculate.sample.prev(data, status, tag.healthy))
	}
	
	cObservedPrev <- measures.acc$cutoffs[which(round(difference,10) == round(min(difference),10))]

	optimal.cutoff <- obtain.optimal.measures(cObservedPrev, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff)

	res
}
