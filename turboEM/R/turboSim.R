## combines two or more objects of class turbo or turbosim into a single turbosim object
## this function should only be called within turboSim(), as it assumes that c("fixptfn", "objfn", "method", "pconstr", "project") remain constant across all objects being combined
combine <- function(..., method.names, keep.pars = FALSE, check = TRUE) {
	## ... = turbo or turbosim objects to be combined
	## keep.pars = logical indicating whether parameter values should be kept
	## check = logical indicating whether the objects should be checked for compatability before being combined
	
	turbo.list <- list(...)
    
    const.vars <- c("fixptfn", "objfn", "method", "pconstr", "project", "control.run", "control.method")
    
    ##tmp <<- turbo.list
    
    ## check that all of the "const.vars" are the same in each object to be combined
    if(check) {
		const.parts <- lapply(turbo.list, function(x) x[const.vars])
	    all.const <- all(sapply(2:length(const.parts), function(x) identical(const.parts[[1]], const.parts[[x]])))
	    if(!all.const) {
	    	stop("Objects should not be combined. Check that the following components are identical across objects: fixptfn, objfn, method, pconstr, project, control.run, and control.method\n")
	    }
	}
    
    combo <- turbo.list[[1]][const.vars]
    
    if(missing(method.names)) {
    	if(!is.null(turbo.list[[1]]$method.names)) {
	    	method.names <- turbo.list[[1]]$method.names
	    } else {
    		method.names <- make.unique(combo$method)
    	}
    }
    combo$method.names <- method.names

    vary.vars <- c("fail", "value.objfn", "itr", "fpeval", "objfeval", "convergence", "errors")

	## for all of the objects that come from a single iteration of the benchmark study, turn each component into a matrix
	## just keep the elapsed time of runtime
    just1fit <- sapply(turbo.list, function(x) is.null(dim(x$fail)))
    for(i in seq_along(turbo.list)) {
    	if(just1fit[i]) {
    		for(j in seq_along(vary.vars)) {
    			turbo.list[[i]][[vary.vars[j]]] <- matrix(turbo.list[[i]][[vary.vars[j]]], nrow=1)
    		}
    		turbo.list[[i]][["runtime"]] <- matrix(turbo.list[[i]][["runtime"]][,"elapsed"], nrow=1)
    	}    		
    }
    vary.vars <- c(vary.vars,"runtime")
    	
	## combine the results across iterations
    for(j in seq_along(vary.vars)) {
    	combo[[vary.vars[j]]] <- do.call("rbind", lapply(turbo.list, function(x) x[[vary.vars[j]]]))
    	colnames(combo[[vary.vars[j]]]) <- method.names
    }
    
    if(keep.pars) {
    	# tmp <- as.vector(lapply(turbo.list, function(x) x[["pars"]]))
    	# combo$pars <- array(unlist(tmp), dim=c(nrow(tmp[[1]]), ncol(tmp[[1]]), length(tmp)), dimnames=list(method=method.names, par=paste("p",1:ncol(tmp[[1]]),sep=""), iter=1:length(tmp)))
    	tmp <- as.vector(lapply(turbo.list, function(x) x[["pars"]]))
    	##niter <- sum(sapply(tmpXXX, function(x) ifelse(is.na(dim(x)[3]), 1, dim(x)[3])))
    	niter <- nrow(combo[[vary.vars[1]]])
    	combo$pars <- array(unlist(tmp), dim=c(nrow(tmp[[1]]), ncol(tmp[[1]]), niter), dimnames=list(method=method.names, par=paste("p",1:ncol(tmp[[1]]),sep=""), iter=1:niter))
    }
    
    # class(combo) <- "turbosim"
    combo
}

turboSim <- function(parmat, fixptfn, objfn, method = c("em", "squarem", "pem", "decme", "qn"), boundary, pconstr = NULL, project = NULL, parallel = FALSE, method.names, keep.pars = FALSE, ..., control.method = replicate(length(method),list()), control.run = list()) {
	combine.nocheck <- function(...) {
		combine(..., method.names=method.names, keep.pars=keep.pars, check=FALSE)
	}
	
	if(parallel) {
		results <- foreach(par = iterators::iter(parmat, by="row"), .packages="turboEM", .combine=combine.nocheck) %dopar% turboem(par=par, fixptfn=fixptfn, objfn=objfn, method=method, boundary=boundary, pconstr=pconstr, project=project, parallel=FALSE, ..., control.method=control.method, control.run=control.run) 
	} else {
		results <- foreach(par = iterators::iter(parmat, by="row"), .packages="turboEM", .combine=combine.nocheck) %do% turboem(par=par, fixptfn=fixptfn, objfn=objfn, method=method, boundary=boundary, pconstr=pconstr, project=project, parallel=FALSE, ..., control.method=control.method, control.run=control.run) 		
	}
	
	if(!missing(method.names)) results$method.names <- method.names
	
	class(results) <- "turbosim"
	results	
}

# turboSim2 <- function(parmat, fixptfn, objfn, method = c("em", "squarem", "pem", "decme", "qn"), boundary, pconstr = NULL, project = NULL, parallel = FALSE, method.names, ..., control.method = replicate(length(method),list()), control.run = list()) {
	# if(parallel) {
		# nworkers <- getDoParWorkers()
		# cuts <- cut(1:nrow(parmat), nworkers, labels=F)
		# results <- foreach(j = 1:nworkers, .packages="turboEM", .combine=combine) %dopar% {
			# foreach(p.start = iter(parmat[cuts==j,], by="row"), .packages="turboEM", .combine=combine) %do% turboem(par=p.start, fixptfn=fixptfn, objfn=objfn, method=method, boundary=boundary, pconstr=pconstr, project=project, parallel=FALSE, ..., control.method=control.method, control.run=control.run)
		# }
	# } else {
		# results <- foreach(p.start = iter(parmat, by="row"), .packages="turboEM", .combine=combine) %do% turboem(par=p.start, fixptfn=fixptfn, objfn=objfn, method=method, boundary=boundary, pconstr=pconstr, project=project, parallel=FALSE, ..., control.method=control.method, control.run=control.run) 		
	# }
	
	# if(!missing(method.names)) results$method.names <- method.names
	
	# class(results) <- "turbosim"
	# results	
# }
