function.MCT <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (control$CFN <= 0 || control$CFP <= 0) {
		stop("You have entered an invalid value for costs. Costs must be positive.", call. = FALSE)
	}
	MCT <- (control$CFN/control$CFP)*pop.prev*(1-measures.acc$Se[,1])+(1-pop.prev)*(1-measures.acc$Sp[,1])
	
	optimal.MCT <- round(min(MCT,na.rm=TRUE),10)	 
	cMCT <- measures.acc$cutoffs[which(round(MCT,10) == optimal.MCT)]
	
	optimal.cutoff <- obtain.optimal.measures(cMCT, measures.acc)
	
	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = MCT, optimal.criterion = optimal.MCT)			
}
