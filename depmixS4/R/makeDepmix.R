

# the main function constructing a mix model with full information, ie all models already in place
# this function is probably not ever called by users ...

makeMix <-
function(response, prior, ...) {
		
	nstates <- length(response)
	nresp <- length(response[[1]])
		
	# count the number of parameters	
	npars <- npar(prior) 
	for(i in 1:nstates) {
		npars <- npars + sum(sapply(response[[i]],npar))
	}
	
	# make appropriate array for response densities
	nt <- nrow(response[[1]][[1]]@y)
	ntimes <- rep(1,nt)
	dens <- array(,c(nt,nresp,nstates))
	
	# compute observation and transition densities
	for(i in 1:nstates) {
		for(j in 1:nresp) {
			dens[,j,i] <- dens(response[[i]][[j]]) # remove this response as an argument from the call to setpars
		}
	}
	
	# compute initial state probabilties
	init <- dens(prior)
	
	new("mix",response=response,prior=prior,
		dens=dens,init=init,nstates=nstates,
		nresp=nresp,ntimes=ntimes,npars=npars)
	
}



# the main function constructing a depmix model with full information, ie all models already in place
# this function is probably not ever called by users

makeDepmix <-
function(response, transition, prior, ntimes=NULL, homogeneous=TRUE, stationary=NULL, ...) {	
	
	if(!(is.null(stationary))) {
			homogeneous <- stationary
			warning("Argument 'stationary' has been replaced by argument 'homogeneous' in 
					version 1.3-0. In future versions argument 'stationary' will most likely be
					used for other purposes.")
	}
	
	nstates <- length(response)
	nresp <- length(response[[1]])
	
	# make appropriate ntimes
	if(is.null(ntimes)) {
		ntimes <- nrow(response[[1]][[1]]@y)
	}
	
	# count the number of parameters	
	npars <- npar(prior) 
	for(i in 1:nstates) {
		npars <- npars + sum(sapply(response[[i]],npar))
	}
	npars <- npars + sum(sapply(transition,npar))
	
	# make appropriate array for transition densities
	nt <- sum(ntimes)
	if(homogeneous) trDens <- array(0,c(1,nstates,nstates))
	else trDens <- array(0,c(nt,nstates,nstates))
	
	# make appropriate array for response densities
	dens <- array(,c(nt,nresp,nstates))
	
	# compute observation and transition densities
	for(i in 1:nstates) {
		for(j in 1:nresp) {
			dens[,j,i] <- dens(response[[i]][[j]]) # remove this response as an argument from the call to setpars
		}
		trDens[,,i] <- dens(transition[[i]])
	}
	
	# compute initial state probabilties
	init <- dens(prior)
	
	# check if dim(init) agrees ntimes
	if(!(dim(init)[1]==length(ntimes))) stop("Argument 'ntimes' does not agree with dimension of prior model.")
	
	new("depmix",response=response,transition=transition,prior=prior,
		dens=dens,trDens=trDens,init=init,homogeneous=homogeneous,
		ntimes=ntimes,nstates=nstates,nresp=nresp,npars=npars)
	
}

