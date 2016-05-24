makeResponseModels <-
function(response,data=NULL,nstates,family,values=NULL,prob=TRUE,...) {
	
	resp <- response
	response <- list()
	nresppars <- 0
		
	# univariate response data
	if(class(resp)=="formula") {
		resp <- list(resp)
		family <- list(family)
	}
	
	if(!(length(resp)==length(family))) stop("Length of 'response' list and 'family' list do not match")
	
	# make response model
	nresp <- length(resp)
	for(i in 1:nstates) {
		response[[i]] <- list()
		for(j in 1:nresp) {
			response[[i]][[j]] <- GLMresponse(resp[[j]],data=data,family=family[[j]])
			nresppars <- nresppars + npar(response[[i]][[j]])
		}
	}
	
	# set the starting values, if any
	if(!is.null(values)) {
		if(!(length(values)==nresppars)) stop(paste("'respstart' has incorrect length, it should be", nresppars, "\n"))
		for(i in 1:nstates) {
			for(j in 1:nresp) {
				bp <- npar(response[[i]][[j]])
				response[[i]][[j]] <- GLMresponse(resp[[j]],data=data,family=family[[j]],pstart=values[1:bp],prob=prob)
				bp <- bp+1
				values <- values[bp:length(values)]
			}
		}
	}
	
	return(response)
}

