function.MaxProdNPVPPV <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	prod <- measures.acc$NPV[,1] * measures.acc$PPV[,1]	 
	cmaxProdNPVPPV <- measures.acc$cutoffs[which(round(prod,10) == round(max(prod,na.rm=TRUE),10))]
	optimal.prod <- max(prod,na.rm=TRUE)

	optimal.cutoff <- obtain.optimal.measures(cmaxProdNPVPPV, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = prod, optimal.criterion = optimal.prod)
	res
}
