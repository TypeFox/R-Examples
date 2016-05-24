#' @title Combining and dividing runjags and MCMC objects
#' @name combine.mcmc
#' @aliases combine.mcmc combine.MCMC combine.jags combine.JAGS divide.jags divide.JAGS
#' @export

#' @description
#' Utility functions for combining separate MCMC or runjags objects into a single object, or the reverse operation

#' @details
#' The combine.mcmc function allows an MCMC object (with 1 or more chains) to be combined with object(s) representing extensions of the same simulation, to produce one MCMC object that contains the continuous combined Markov chains.  Alternatively, a single MCMC list object can be converted into a single chain by combining all chains sequentially.  An object of class \code{\link{runjags-class}} can also be used, in which case the MCMC objects will be extracted from this.  The combine.jags function does a similar operation, but returning the entire runjags object as a single object that can be extended using \code{\link{extend.jags}}.  The divide.jags extracts one or more chains from a given runjags object.

#' @keywords methods

#' @return
#' For combine.mcmc:  an MCMC object if collapse.chains=TRUE, or an mcmc.list object if collapse.chains=FALSE
#'
#' For combine.jags and divide.jags:  a \code{\link{runjags-class}} object

#' @seealso
#' \code{\link{run.jags}} and \code{\link{runjags-class}}

#' @param mcmc.objects a list of MCMC or runjags objects, all with the same number of chains and matching variable names, or a single MCMC object/list or runjags object.  No default.

#' @param thin an integer to use to thin the (final) MCMC object by, in addition to any thinning already applied to the objects before being passed to combine.mcmc.  Ignored if return.samples is specified (!is.na).  Default 1 (no additional thinning is performed).

#' @param return.samples the number of samples to return after thinning.  The chains will be thinned to as close to this minimum value as possible, and any excess iterations discarded.  Supersedes thin if  both are specified.  Ignored if niter(mcmc.objects) < return.samples.   Default NA.

#' @param collapse.chains option to combine all MCMC chains into a single MCMC chain with more iterations.  Can be used for combining chains prior to calculating results in order to reduce the Monte Carlo error of estimates.  Default TRUE if a single mcmc.object is provided, or FALSE otherwise.

#' @param vars an optional character vector of variable names to extract.  If supplied, only variable names in the object supplied with a partial match to anything in 'vars' will be used.  Note that regular expressions are not allowed, but the caret (^) token can be used to specify the match at the start of a variable name, and a quoted vars will be matched exactly.  Default NA meaning all variables available are returned.

#' @param add.mutate should any mutate function associated with the runjags objects be run to collect the additional variables before returning MCMC chains?

#' @param runjags.objects a list of runjags class objects to be combined

#' @param runjags.object a single runjags class object to be divided

#' @param which.chains the chains to extract from the runjags object

#' @param summarise option to add a new set of summary statistics to the newly created runjags object

#' @param ... other arguments to be passed to \code{\link{add.summary}}
NULL


#' @rdname combine.mcmc
combine.mcmc <- function(mcmc.objects=list(), thin=1, return.samples=NA, collapse.chains=if(length(mcmc.objects)==1) TRUE else FALSE, vars=NA, add.mutate=TRUE){

	if(class(mcmc.objects)!="list"){
		if(inherits(mcmc.objects,c("mcmc.list", "mcmc", "runjags"))){
			mcmc.objects <- list(mcmc.objects)
		}else{
			stop("Data must be provided as a list of or single mcmc object(s), or a list of or single mcmc.list(s) (for multiple chains)")
		}
	}
	
	if(!class(thin)%in%c('numeric','integer')) stop('Invalid value provided for "thin"')
	
	mcmc.objects <- lapply(mcmc.objects, function(x){
		if(class(x)=="runjags") return(as.mcmc.list(x, add.mutate=add.mutate)) else return(x)   # using as.mcmc.list adds the mutated function
	})
		
	returnmcmcs <- all(sapply(mcmc.objects,class)=="mcmc")
	if(returnmcmcs && collapse.chains){
		warning("Can't collapse chains of single-chain mcmc objects")
		collapse.chains <- FALSE
	}
	
	mcmc.objects <- lapply(mcmc.objects, function(x){
		if(class(x)=="mcmc") return(as.mcmc.list(x)) else return(x)
	})
		
	no.objects <- length(mcmc.objects)
	
	if(length(mcmc.objects)==0)
		stop("The list provided cannot be empty")
	if(any(sapply(mcmc.objects, niter)==0))
		stop("One or more mcmc objects has 0 iterations")
	
	n.chains <- integer(length=no.objects)
	n.params = rowlengths <- vector("list", length=no.objects)
	
	# mcmc.objects is guaranteed to be a list of mcmc.lists:
	vnames <- lapply(mcmc.objects, varnames)
	# Check that all varnames are the same
	if(length(mcmc.objects)>1){
		allsame <- sapply(vnames, function(x) return(all(x==vnames[[1]])))
		if(!all(allsame)) stop("Non matching variable names for supplied MCMC objects")
	}
		
	# Check now that we have at least 1 matching variable name:
	fornames <- mcmc.objects[[1]][[1]]
	selected <- matchvars(checkvalidmonitorname(vars),  varnames(fornames))
				
	if(is.null(dimnames(fornames)) || is.null(dimnames(fornames)[[1]])) iterstart <- 1 else iterstart <- as.numeric(dimnames(fornames)[[1]])[1]
	iterthin <- thin(fornames)
	
	if(no.objects > 1){
		
		for(i in 1:no.objects){
		
			if(class(mcmc.objects[[i]])=="mcmc.list"){
				n.chains[i] <- length(mcmc.objects[[i]])
				returnlist <- TRUE
			}else{
				if(class(mcmc.objects[[i]])=="mcmc"){
					n.chains[i] <- 1
					mcmc.objects[[i]] <- mcmc.list(mcmc.objects[[i]])
					returnlist <- FALSE
				}else{
					stop("Data must be provided as a list of mcmc objects, or a list of mcmc.lists (for multiple chains).")	
				}
			}
		
			n.params[[i]] <- integer(length=n.chains[i])
			rowlengths[[i]] <- integer(length=n.chains[i])
		
			for(j in 1:n.chains[i]){
				n.params[[i]][j] <- nvar(mcmc.objects[[i]][[j]])
				rowlengths[[i]][j] <- niter(mcmc.objects[[i]][[j]])
			}
		
			if(!all(rowlengths[[i]] == rowlengths[[i]][1])) stop(paste("The chain lengths were not equal for object ", i, ".", sep=""))
			rowlengths[[i]] <- rowlengths[[i]][1]
			
		}
	
		rowlengths <- unlist(rowlengths)

		paramsequal <- all(unlist(lapply(n.params, function(x) if(all(x==n.params[[1]][1])) return(TRUE) else return(FALSE))))
	
		if(!(all(n.chains==n.chains[1]))) ("There was an unequal number of chains between mcmc objects")
		if(!paramsequal) stop("There was an unequal number of monitored variables (columns) between chains / mcmc objects")
	
		n.chains <- n.chains[1]
		n.params <- n.params[[1]][1]
	
		newobjects <- vector("list", length=n.chains)
		
		for(i in 1:n.chains){
		
			newobjects[[i]] <- matrix(NA, nrow=0, ncol=n.params, dimnames=list(NULL, dimnames(mcmc.objects[[1]][[1]])[[2]]))
		
			for(j in 1:no.objects){
				
				if(is.null(dim(mcmc.objects[[j]][[i]]))) dim(mcmc.objects[[j]][[i]]) <- c(niter(mcmc.objects[[j]][[i]]), 1)
				newobjects[[i]] <- rbind(newobjects[[i]], mcmc.objects[[j]][[i]])
			
			}
			
			rowlengths <- nrow(newobjects[[i]])
			newiternames <- (((iterstart-1)/iterthin):((((iterstart-1)/iterthin)+rowlengths)-1) * iterthin) + 1
			
			dimnames(newobjects[[i]])[[1]] <-  newiternames
			
			newobjects[[i]] <- mcmc(newobjects[[i]], start=iterstart, thin=iterthin)

		}
		
		if(returnlist) newobjects <- as.mcmc.list(newobjects) else newobjects <- newobjects[[1]]
		
	}else{
		
		newobjects <- mcmc.objects[[1]]
	}
	
	rowlengths <- niter(newobjects)

	# For combine.mcmc part later, more efficient to thin out first if we are doing that anyway (recursive call will thin further as required):
	startretsamples <- return.samples
	startthin <- thin	
		
	if(!is.na(return.samples)){
		if(return.samples > rowlengths){
			thin <- 1
			if(return.samples!=Inf) warning('Specified return.samples was longer than the chains provided - returning shorter MCMC object length')
		}else{
			thin <- rowlengths / return.samples
		}
	}else{
		return.samples <- Inf
	}

	currentthin <- thin(newobjects)
	thin <- floor(thin)*currentthin
	endat <- (stats::start(newobjects)+(thin*return.samples))-1

	suppressWarnings(newobjects <- window(newobjects, thin=thin, end=endat))

	# And only spit back the variables requested (checked at the start):
	thevarnames <- dimnames(newobjects[[1]])
	newobjects <- newobjects[,selected,drop=FALSE]
	for(i in 1:length(newobjects)){
		dimnames(newobjects[[i]]) <- list(thevarnames[[1]], thevarnames[[2]][selected])
	}
	
	if(is.null(dimnames(newobjects[[1]])) && is.null(dimnames(newobjects[[1]])[[1]])){
		warning("NULL iteration names produced")
	}
	
	# If collapse chains recursive call (after thinning to approx double what we will need, for efficiency):
	if(collapse.chains){
		class(newobjects) <- "list"
		newobjects <- combine.mcmc(newobjects, collapse.chains=FALSE, return.samples=startretsamples, thin=1, vars=NA)	
	}else{
		if(returnmcmcs) newobjects <- newobjects[[1]]			
	}
		
	return(newobjects)
	
}

combine.MCMC <- combine.mcmc


#' @rdname combine.mcmc
combine.jags <- function(runjags.objects=list(), summarise=TRUE, ...){
	
	if(class(runjags.objects)!='list' || length(runjags.objects)==0)
		stop('A list of runjags objects must be supplied')
	
	if(length(runjags.objects)==1)
		stop('Only 1 runjags object provided - cannot combine!')
	
	runjags.objects <- lapply(runjags.objects, checkvalidrunjagsobject)
	
	if(!all(runjags.objects[[1]]$model == sapply(runjags.objects, function(x) return(x$model))))
		stop('The list of runjags objects supplied must have identical model formulations')
	if(!all(runjags.objects[[1]]$data == sapply(runjags.objects, function(x) return(x$data))))
		stop('The list of runjags objects supplied must have identical data')	
	if(any(sapply(runjags.objects, function(x) return(any(x$modules!='')))) && all(paste(runjags.objects[[1]]$modules,collapse='\n') == sapply(runjags.objects, function(x) return(paste(x$modules,collapse='\n')))))
		stop('The list of runjags objects supplied must have identical modules')	
	if(any(sapply(runjags.objects, function(x) return(any(x$factories!='')))) && all(paste(runjags.objects[[1]]$factories,collapse='\n') == sapply(runjags.objects, function(x) return(paste(x$factories,collapse='\n')))))
		stop('The list of runjags objects supplied must have identical factories')	
	
	summaryargs <- list(...)
	summaryargs <- getargs('add.summary', summaryargs, otherfnames='combine.jags')
	
	monitors.same <- FALSE
	if(all(runjags.objects[[1]]$monitor == sapply(runjags.objects, function(x) return(x$monitor))))
		monitors.same <- TRUE

	chains.same <- FALSE
	if(all(paste(runjags.objects[[1]]$end.state,collapse='\n') == sapply(runjags.objects, function(x) return(paste(x$end.state,collapse='\n')))))
		chains.same <- TRUE
	
	if(chains.same)# && monitors.same)
		stop('All runjags objects supplied appear to be identical')
		
#	if(!monitors.same && !chains.same)
#		stop('The runjags objects to combined must either have identical monitored variables (to combine chains), or identical end states (to combine monitors that have been read in separately)')
	
	if(!monitors.same)
		stop('The runjags objects to combined must have identical monitored variables')
	
	# Use the maximum time taken:
	ntt <- sapply(runjags.objects, function(x) return(max(as.double(x$timetaken, units='secs'))))
	timetaken <- runjags.objects[[order(ntt,decreasing=TRUE)[1]]]$timetaken
	
	if(monitors.same){
		# Combine chains and end.state
		newobj <- runjags.objects[[1]]
		
		newobj$mcmc <- as.mcmc.list(unlist(lapply(runjags.objects, function(x) return(x$mcmc)), recursive=FALSE))   # This will NOT add the mutate function
		newobj$end.state <- unlist(lapply(runjags.objects, function(x) return(x$end.state)))
		class(newobj$end.state) <- 'runjagsinits'
	}
	
	# Maybe implement at some point:
#	if(chains.same){		
#		newobj <- runjags.objects[[1]]
#		allvarnames <- lapply(runjags.objects, function(x) return(varnames(as.mcmc.list(x))))
#		# combine monitors and work out what should be noread.monitor
#	}

	newobj <- makerunjagsobject(newobj, summarise, summaryargs, newobj$burnin, newobj$sample, newobj$thin, newobj$model, newobj$data, newobj$monitor, newobj$noread.monitor, newobj$modules, newobj$factories, response=newobj$response, residual=newobj$residual, fitted=newobj$fitted, newobj$method, newobj$method.options, timetaken, silent=FALSE)
		
	return(newobj)
}

combine.JAGS <- combine.jags

#' @rdname combine.mcmc
divide.jags <- function(runjags.object, which.chains=1:nchain(as.mcmc.list(runjags.object)), summarise=TRUE, ...){
	
	# Add vars as an option?
	
	runjags.object <- checkvalidrunjagsobject(runjags.object)
	
	if(is.na(which.chains) || any(which.chains > nchain(as.mcmc.list(runjags.object))))
		stop('Invalid value provided for the which.chain argument')
	
	summaryargs <- list(...)
	summaryargs <- getargs('add.summary', summaryargs, otherfnames='divide.jags')
	
	newobj <- runjags.object
	
	newobj$mcmc <- runjags.object$mcmc[which.chains]
	newobj$end.state <- runjags.object$end.state[which.chains]
	
	newobj <- makerunjagsobject(newobj, summarise, summaryargs, newobj$burnin, newobj$sample, newobj$thin, newobj$model, newobj$data, newobj$monitor, newobj$noread.monitor, newobj$modules, newobj$factories, response=newobj$response, residual=newobj$residual, fitted=newobj$fitted, newobj$method, newobj$method.options, newobj$timetaken, silent=FALSE)
		
	return(newobj)
}

divide.JAGS <- divide.jags

