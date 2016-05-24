makePriorModel <-
function(nstates,ncases,formula=~1,data=NULL,values=NULL, ...) {
	
	# these arguments need to be added at some point FIX ME
	base=1
	# prob=TRUE
	
	# initial probabilities model, depending on covariates init(=~1 by default)
	if(formula==~1) {
		initModel <- transInit(~1,data=data.frame(rep(1,ncases)),nstates=nstates,family=multinomial(link="identity"),pstart=values)
	} else {
		if(is.null(data)) {
			stop("Argument 'initdata' missing with non-trivial prior model")
		} else {
			initModel <- transInit(formula,data=data,nstates=nstates,family=multinomial(),pstart=values)
		}	
	}
	
	return(initModel)
}

