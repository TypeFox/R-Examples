function.MaxNPVPPV <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	m <- vector()
	for(i in 1:length(measures.acc$cutoffs)) {
		if (measures.acc$PPV[i,1] <= measures.acc$NPV[i,1]) {
			m[i] <- measures.acc$PPV[i,1]
		} else  {
			m[i] <- measures.acc$NPV[i,1]
		}
	}
	M <- max(m,na.rm=TRUE)   

	optimal.index <- which(round(m,10) == round(M,10))
	cMaxNPVPPV <- measures.acc$cutoffs[optimal.index]

	optimal.cutoff <- obtain.optimal.measures(cMaxNPVPPV, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = m, optimal.criterion = M)
	res
}
