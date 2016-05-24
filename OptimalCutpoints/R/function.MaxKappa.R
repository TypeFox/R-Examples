function.MaxKappa <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)  
	if (is.logical(control$weighted.Kappa) == FALSE) {
		stop("'weighted.Kappa' must be a logical-type argument.", call. = FALSE)
	}
	if (control$weighted.Kappa == FALSE) {
                costs.rate <- 0.5		
	}
	if (control$weighted.Kappa == TRUE) {
		if (control$CFN <= 0 || control$CFP <= 0) {
			stop("You have entered an invalid value for costs. Costs must be positive.", call. = FALSE)
		}		
		costs.rate <- control$CFN/(control$CFP+control$CFN)
	}		
	
	Kappa <-(pop.prev*(1-pop.prev)*(measures.acc$Se[,1]+measures.acc$Sp[,1]-1))/(pop.prev*(pop.prev*(1-measures.acc$Se[,1])+(1-pop.prev)*measures.acc$Sp[,1])*costs.rate+(1-pop.prev)*(pop.prev*measures.acc$Se[,1]+
		(1-pop.prev)*(1-measures.acc$Sp[,1]))*(1-costs.rate))
	
	cMaxKappa <- measures.acc$cutoffs[which(round(Kappa,10) == round(max(Kappa, na.rm = TRUE),10))]
	optimal.Kappa <- max(Kappa, na.rm = TRUE)
	optimal.cutoff <- obtain.optimal.measures(cMaxKappa, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = Kappa, optimal.criterion = optimal.Kappa)
	res
}