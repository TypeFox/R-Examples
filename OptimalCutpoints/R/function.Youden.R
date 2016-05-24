function.Youden <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	if (is.logical(control$generalized.Youden) == FALSE) {
		stop("'generalized.Youden' must be a logical-type argument.", call. = FALSE)		
	}
	if (is.logical(control$costs.benefits.Youden) == FALSE) {
		stop("'costs.benefits.Youden' must be a logical-type argument.", call. = FALSE)
	}		 
	if (control$generalized.Youden == FALSE) {				
		expression.Youden <- measures.acc$Se[,1] + measures.acc$Sp[,1]-1		 
	}
	if (control$generalized.Youden == TRUE) {
		if (control$CFN <= 0 || control$CFP <= 0) {
			stop("You have entered an invalid value for costs. Costs must be positive.", call. = FALSE)
		}
		r <- ((1-pop.prev)/pop.prev)*(control$CFP/control$CFN)
		expression.Youden <- measures.acc$Se[,1]+r*measures.acc$Sp[,1]-1
	}

	if (control$costs.benefits.Youden == FALSE) {
			cYouden <- measures.acc$cutoffs[which(round(expression.Youden,10) == round(max(expression.Youden, na.rm=TRUE),10))]			
	}
	
	if (control$costs.benefits.Youden == TRUE & control$generalized.Youden == FALSE) {
		control$costs.ratio <- 1
		pop.prev <- 0.5
		cYouden <- function.CB(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control, pop.prev, ci.fit, conf.level, measures.acc)$optimal.cutoff$cutoff
	}
	
	if (control$costs.benefits.Youden == TRUE & control$generalized.Youden == TRUE) {
		control$costs.ratio = control$CFP/control$CFN

		cYouden <- function.CB(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control, pop.prev, ci.fit, conf.level, measures.acc)$optimal.cutoff$cutoff
	}
		
	optimal.cutoff <- obtain.optimal.measures(cYouden, measures.acc)
	
	if (control$generalized.Youden == FALSE) {
		Youden <- unique(round(optimal.cutoff$Se[,1]+optimal.cutoff$Sp[,1]-1, 10))
	}
	
	if (control$generalized.Youden == TRUE) {
		Youden <- unique(round(optimal.cutoff$Se[,1]+r*optimal.cutoff$Sp[,1]-1, 10))
	}
	
	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = expression.Youden, optimal.criterion = Youden)		
	res
}
