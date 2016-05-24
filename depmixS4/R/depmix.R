#
# Ingmar Visser, 11-6-2008
#

#
# Main function to construct mix models
#

#
# UNIVARIATE AND MULTIVARIATE MIXTURE OF GLM'S
#


setGeneric("mix", function(response, data = NULL, 
    nstates, family = gaussian(), prior = ~1, initdata = NULL, 
    respstart = NULL, instart = NULL, ...) standardGeneric("mix"))


setMethod("mix", signature(response = "ANY"), function(response, 
    data = NULL, nstates, family = gaussian(), prior = ~1, initdata = NULL, 
    respstart = NULL, instart = NULL, ...) {
    
	#if(!is.null(data) & any(is.na(data))) stop("'depmixS4' does not currently handle missing data.")
	
    # make response models
    response <- makeResponseModels(response = response, data = data, 
        nstates = nstates, family = family, values = respstart)
    
    # FIX ME: this only works if data are actually provided ... 
	# (maybe make this obligatory ...)
    ntimes <- rep(1, nrow(data))
    
    # make prior model
	if(is.null(initdata)) initdata=data
    prior <- makePriorModel(nstates = nstates, ncases = length(ntimes), 
        formula = prior, data = initdata, values = instart)
    
    # call main depmix with all these models, ntimes and homogeneous
    model <- makeMix(response = response, prior = prior)
        
    return(model)
})

#
# Ingmar Visser, 23-3-2008
#

#
# Main function to construct depmix models
#

#
# UNIVARIATE AND MULTIVARIATE MARKOV MIXTURE OF GLM'S
#

setMethod("depmix", signature(response = "ANY"), function(response, 
		data = NULL, nstates, transition = ~1, family = gaussian(), 
		prior = ~1, initdata = NULL, respstart = NULL, trstart = NULL, 
		instart = NULL, ntimes = NULL, ...) {
    
	if(is.null(data)) {
		if(is.null(ntimes)) stop("'ntimes' must be provided if not in the data")
	} else {
		#if(any(is.na(data))) stop("'depmixS4' does not currently handle missing data.")
		if(is.null(attr(data, "ntimes"))) {
			if (is.null(ntimes)) ntimes <- nrow(data)
		} else {
			ntimes <- attr(data, "ntimes")
		}
		if (sum(ntimes) != nrow(data)) stop("'ntimes' and data do not match")
	}
	
	if(nstates==1&transition!=~1) {stop("1-state model can not have transition covariate")}
	
	# make response models
    response <- makeResponseModels(response = response, data = data, 
        nstates = nstates, family = family, values = respstart)
    
    # make transition models
    homogeneous = FALSE
    if (transition == ~1) 
        homogeneous = TRUE
    transition <- makeTransModels(nstates = nstates, formula = transition, 
        data = data, homogeneous = homogeneous, values = trstart)
    
    # make prior model
    prior <- makePriorModel(nstates = nstates, ncases = length(ntimes), 
        formula = prior, data = initdata, values = instart)
    
    # call main depmix with all these models, ntimes and homogeneous
    model <- makeDepmix(response = response, transition = transition, 
        prior = prior, ntimes = ntimes, homogeneous = homogeneous)
    
    # deal with starting values here!!!!!!
    
    return(model)
})
