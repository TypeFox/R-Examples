uncertainty <- 
function(aep, uc.values, uc.names, prob=seq(5,95,5), digits=c(0,0), print=TRUE) {
###	uncertainty assessment of annual energy production
	
	if(missing(aep)) stop("AEP object 'aep' is mandatory")
	if(missing(uc.values)) stop("Uncertainty values 'uc.values' is mandatory")
	if(missing(uc.names)) stop("Uncertainty names 'uc.names' is mandatory")
	if(class(aep)!="aep") stop(substitute(aep), " is no aep object")
	if(!is.numeric(uc.values)) stop("'uc.values' must be numeric")
	for(i in 1:length(uc.values)) if(uc.values[i]<0) stop("Only positive 'uc.values' allowed")
	if(!is.null(uc.names)) if(length(uc.names)!=length(uc.values) && length(uc.names)!=length(uc.values)+1) stop("'uc.names' and 'uc.values' must be vectors of the same length")
	if(length(digits)!=2) { 
		if(length(digits)<2) digits <- rep(digits, 2)
		if(length(digits)>2) digits <- digits[1:2]
		warning("'digits' shall be a vector of two values (for uncertainty and AEP values)", call.=FALSE)
	}
	
	p50 <- tail(aep$aep$total, 1)
	uc.tot <- sqrt(sum(uc.values^2))
	prob.ex <- data.frame(cbind(prob, qnorm(rev(prob / 100), p50, uc.tot / 100 * p50)))
	names(prob.ex) <- c("probability", "aep")
	prob.ex$aep <- round(prob.ex$aep, digits[2])
	attr(prob.ex$probability, "unit") <- "%"
	attr(prob.ex$aep, "unit") <- attr(aep$aep$total, "unit")
	attr(prob.ex$aep, "P50") <- p50
	
	uc <- data.frame(c(uc.values, uc.tot))
	names(uc) <- "uncertainty"
	if(is.null(uc.names)) row.names(uc)[length(uc$uncertainty)] <- "total" 
	else {
		if(length(uc.names)==length(uc.values))	row.names(uc) <- c(uc.names, "total")
		else row.names(uc) <- uc.names
	}
	uc <- round(uc, digits[1])
	attr(uc, "unit") <- "%"
	
	uncertainty <- list(uncertainty.meth=uc, prob.exceedance=prob.ex)
	
	attr(uncertainty, "call") <- list(func="uncertainty", aep=deparse(substitute(aep)), uc.values=uc.values, uc.names=uc.names, prob=prob, digits=digits, print=print)
	class(uncertainty) <- "uncertainty"
	
	if(print) print(uncertainty)
	invisible(uncertainty)	
}
