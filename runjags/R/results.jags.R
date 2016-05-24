#' @title Importing of saved JAGS simulations with partial error recovery
#' @name results.jags
#' @aliases results.jags results.JAGS
#' @export

#' @description
#' Imports a completed JAGS simulation from a folder created by \code{\link{run.jags}} using the background or bgparallel methods, or any other method where the keep.jags.files=TRUE option was used. Partial recovery simulations is possible for parallel methods where one or more simulation failed to complete. Additional chain thinning and parameter import selection is also supported.

#' @keywords models

#' @return
#' An object of class 'runjags' (see \code{\link{runjags-class}}).

#' @seealso
#' \code{\link{runjags-class}} for details of available methods for the returned object, \code{\link{run.jags}} for details of how to start simulations, and \code{\link{runjags.options}} for user options regarding warning messages etc.

#' @examples

#' \dontshow{
#' runjags.options(new.windows=FALSE)
#' }
#' # Run a model using parallel background JAGS calls:
#' 
#' # Simulate the data:
#' N <- 100
#' X <- 1:N
#' Y <- rnorm(N, 2*X + 10, 1)
#' # Initial values for 2 chains:
#' m <- list(-10, 10)
#' c <- list(-10, 10)
#' precision <- list(0.01, 10)
#'
#' # Model in the JAGS format
#' model <- "model { 
#' for(i in 1 : N){ 
#' 	Y[i] ~ dnorm(true.y[i], precision);
#' 	true.y[i] <- (m * X[i]) + c
#' } 
#' m ~ dunif(-1000,1000)
#' c ~ dunif(-1000,1000) 
#' precision ~ dexp(1)
#' #data# X, Y, N
#' #monitor# m, c, precision
#' #inits# m, c, precision
#' }"
#' 
#' \dontrun{
#' # Run the model and produce plots 
#' fileinfo <- run.jags(model=model, n.chains=2, method="bgparallel")
#' # Wait for the simulations to complete:
#' Sys.sleep(10)
#' # Import only variable m from the first chain:
#' results <- results.jags(fileinfo, read.monitor='m', recover.chains=1)
#' # Look at the summary statistics:
#' print(results)
#' }

#' @param foldername the absolute or relative path to the folder containing the JAGS simulation to be imported.  May also be the return value of a call to \code{\link{run.jags}} with method = 'background' or method = 'bgparallel', which will avoid having to re-load some information from disk and therefore may be slightly faster.

#' @param echo option to display the output of the simulations to screen.  If the simulations have not finished, the progress-to-date will be displayed.

#' @param combine a logical flag indicating if results from the new JAGS run should be combined with the previous chains.  Default value respects the setting chosen during the initial run.jags function call, changing the option to TRUE or FALSE overrides the original setting.

#' @param summarise should summary statistics be automatically calculated for the output chains?  Default value respects the setting chosen during the initial run.jags function call, changing the option to TRUE or FALSE overrides the original setting.

#' @param keep.jags.files option to keep the folder with files needed to call JAGS, rather than deleting it after importing.  Default value respects the setting chosen during the initial run.jags function call, changing the option to TRUE or FALSE overrides the original setting.  See also the \code{\link{cleanup.jags}} function.

#' @param read.monitor an optional character vector of variables to import from the simulation folder.  This may be useful for models with large numbers of variables that would otherwise not be able to be loaded into R.  Default value loads all variables given by the monitor argument (but NOT the noread.monitor argument) to the original run.jags call.

#' @param return.samples option to thin the final MCMC chain(s) before calculating summary statistics and returning the chains.  Note that this option does NOT currently carry out thinning in JAGS, therefore R must have enough available memory to hold the chains BEFORE thinning (for very large chains, it may be necessary to specify a subset of the variables at a time using read.monitor='...' and keep.jags.files=TRUE).  Default value returns all available iterations.

#' @param recover.chains option to try to recover successful simulations if some simulations failed (this is only relevant for parallel methods with more than 1 simulation).  A value of TRUE returns only successful simulations, FALSE will cause an error if any simulation has failed.  A numeric vector of specific chain(s) to be read is also permitted, but an error will be returned if any of the simulations containing these chains was unsuccessful.  The default version reads the option set in \code{\link{runjags.options}}.

#' @param ... additional summary parameters to be passed to \code{\link{add.summary}}
NULL


#' @rdname results.jags
results.jags <- function(foldername, echo=NA, combine=NA, summarise=NA, keep.jags.files=NA, read.monitor=NA, return.samples=NA, recover.chains=NA, ...){
	
	starttime <- Sys.time()
	
	# I switched the argument names:
	sub.chains <- recover.chains    # either FALSE (stop if any crashed), TRUE (try to recover), NA (use runjags.getOption), or a numeric vector of which chains to read in
	sub.samples <- return.samples	# either NA (use sample size) or a number of iterations to thin to
	if(is.na(sub.samples)) sub.samples <- FALSE
		
	background.runjags.object <- foldername
	if(is.character(background.runjags.object) && file.exists(file.path(background.runjags.object, 'jagsinfo.Rsave'))){
		# The folder may have been moved, so save the new path (make it an absolute path):
		newpath <- normalizePath(background.runjags.object, winslash='/')
		# A copy of the original background.runjags.object will be loaded from here:
		load(file.path(background.runjags.object, 'jagsinfo.Rsave'))
		# If it's been moved update the new path:
		if(newpath!=background.runjags.object$directory){
			background.runjags.object$directory <- newpath
			# Give a warning or display of some kind here?
		}
	}
	if(!inherits(background.runjags.object, "runjagsbginfo")){
		if(is.character(background.runjags.object)){
			stop(paste("No valid runjags simulation found at the path provided: '", background.runjags.object, "'", sep=''))
		}else{
			stop("An object produced by a background runjags method (or a path to the JAGS folder to be imported) must be supplied (see the manual page for more details)")
		}
	}
	if(background.runjags.object$method=="xgrid"){
		warning("Re-routing function call to xgrid.results.jags to retrieve xgrid job")
		return(xgrid.results.jags(background.runjags.object))
	}
	background.runjags.object <- checkvalidrunjagsobject(background.runjags.object)
	
	if(!file.exists(background.runjags.object$directory)){
		stop(paste("The JAGS files for the model supplied were not found at the saved directory '", background.runjags.object$directory, "' - try calling results.jags with the path to the new location of the simulation folder", sep=""))
	}
	
	if(!identical(summarise, NA)){
		background.runjags.object$summarise <- summarise
	}
	
	
	if(!identical(combine, NA)){
		if(combine && !background.runjags.object$combine) stop('Can only combine chains if the original function call was to extend.jags with combine=TRUE')
		background.runjags.object$combine <- combine
	}
	combine <- background.runjags.object$combine
	# This was done in extend.jags before, but now we may override the original request for combine so do it here instead:
	if(background.runjags.object$combine) background.runjags.object$burnin <- 0 else background.runjags.object$burnin <- background.runjags.object$burnin + background.runjags.object$adapt
	
	if(identical(sub.samples, TRUE) || identical(sub.samples, NA)) stop('Value for "sub.samples" must be either FALSE or the number of iterations to return')
	
	if(identical(sub.chains, NA)){
		if(!background.runjags.object$combine && runjags.getOption("partial.import")) sub.chains <- TRUE else sub.chains <- FALSE
	}

	if(!identical(sub.chains, FALSE) && background.runjags.object$combine) stop('Cannot return partially completed jobs when combining MCMC objects - try again with combine=FALSE')
	
	
	# Make name easier to type:
	runjags.object <- background.runjags.object
	
	if(!background.runjags.object$summarise && !identical(list(...), list()) && runjags.getOption('summary.warning'))
		warning('Options to add.summary ignored when summarise=FALSE')
	summaryargs <- getsummaryargs(runjags.object$summary.pars, 'results.jags', ...)		
	
	if(!is.na(keep.jags.files)) runjags.object$keep.jags.files <- keep.jags.files
	if(!identical(read.monitor, NA)){
		read.monitor <- checkvalidmonitorname(read.monitor)
		# Switch things in monitor not in read.monitor to noread.monitor:
		tomove <- runjags.object$monitor[! runjags.object$monitor %in% read.monitor]
		# Preserve these (not really monitors anyway):
		tokeep <- runjags.object$monitor[runjags.object$monitor %in% c('pd','full.pd','popt','dic')]
		runjags.object$monitor <- c(read.monitor, tokeep)
		runjags.object$noread.monitor <- c(runjags.object$noread.monitor, tomove)
		# If read.monitor is a subset of an array, the whole array will be moved to noread and the subset kept in monitor (but with indices expanded)
	}
	
	# Call runjags.readin and then deal with copying files etc:
	
	if(length(runjags.object$noread.monitor)>0){
		read.monitor <- runjags.object$monitor[!runjags.object$monitor%in%c('pd','full.pd','popt','dic')]
	}else{
		read.monitor <- NA
	} 
	
	allok <- FALSE
	# Make sure nothing gets deleted if it crashes anywhere:
	reallydelete <- runjags.object$keep.jags.files
	runjags.object$keep.jags.files <- TRUE
	
	# This will be called if runjags.readin fails (allok=FALSE) or at the end of the function (allok=TRUE):
	on.exit({
		
		new.directory <- runjags.object$directory
		
		if(runjags.object$keep.jags.files){
			if(!allok)
			  swcat('Note: Either one or more simulation(s) failed, or there was an error in processing the results.  You may be able to retrieve any successful simulations using:\nresults.jags("', new.directory, '", recover.chains=TRUE)\nSee the help file for that function for possible options.\n', sep='')
			
			# Don't add to the delete on exit list:
#			if(!new.directory %in% runjagsprivate$simfolders)
#				runjagsprivate$simfolders <- c(runjagsprivate$simfolders, new.directory)
		}else{
			if(!allok && runjags.getOption('keep.crashed.files') && new.directory!="Directory not writable"){
			  swcat('Note: Either one or more simulation(s) failed, or there was an error in processing the results.  You may be able to retrieve any successful simulations using:\nresults.jags("', new.directory, '", recover.chains=TRUE)\nSee the help file for that function for possible options.\n', sep='')
			  			  
				swcat('To remove failed simulation folders use cleanup.jags() - this will be run automatically when the runjags package is unloaded\n')
				runjags.object$keep.jags.files <- TRUE
			}
		}
		if(!allok){
			runjagsprivate$failedsimfolders <- c(runjagsprivate$failedsimfolders, new.directory)
			# Remove the failed simulation from the list of sims:
			runjagsprivate$simfolders <- runjagsprivate$simfolders[runjagsprivate$simfolders!=new.directory]			
		}

		# If we either want to delete files or we have copied from the temp dir, delete the sim folder:
		if(!runjags.object$keep.jags.files) unlink(runjags.object$directory, recursive = TRUE)
					
	})
	
	newoutput <- runjags.readin(directory=runjags.object$directory, silent.jags=runjags.object$silent.jags, target.adapt=runjags.object$adapt, target.burnin=(runjags.object$burnin-runjags.object$adapt), target.iters=runjags.object$sample, n.chains=length(runjags.object$inits), monitor=runjags.object$monitor, method=runjags.object$method, method.options=runjags.object$method.options, suspended=TRUE, showoutput=echo, read.monitor=read.monitor, sub.samples=sub.samples, sub.chains=sub.chains, force.read=length(runjags.object$noread.monitor)>0)
	allok <- TRUE
	
	if('errormessage'%in%names(newoutput))
		stop(newoutput$errormessage,call.=FALSE)	
	
	# unfinished will be NULL or TRUE:
	if(identical(newoutput$unfinished,TRUE)){		
		if(length(newoutput$finished.chains==1)){
			errmsg <- paste('The simulation has not finished', sep='')
		}else{
			errmsg <- paste('Simulation(s) ', paste(which(!newoutput$finished.chains), collapse=', '), ' have not finished', sep='')
		}
		runjags.object$keep.jags.files <- TRUE
		if(identical(echo, FALSE)) errmsg <- paste(errmsg, ' - call the same function again with echo=TRUE to see the current progress', sep='')
		stop(errmsg)
	}
	
	simtimetaken <- difftime(newoutput$end.time, runjags.object$startedon)
	
	end.state <- newoutput$end.state
	class(end.state) <- 'runjagsinits'
	
	# Just for iteration names:
	burnin <- runjags.object$oldburnin+(runjags.object$oldthin*runjags.object$oldsample)+runjags.object$burnin
	if(identical(sub.samples, FALSE)){
		# Sample is always not including thinning:
		sample <- runjags.object$sample
		thin <- runjags.object$thin
	}else{
		sample <- min(sub.samples, runjags.object$sample)
		
		# This means that if we extend, the new thin is used - not ideal:
		# thin <- max(1,floor(runjags.object$sample/sub.samples))*runjags.object$thin
		
		# So thin will be wrong for this sample - documented
		thin <- runjags.object$thin
		runjags.object$sample <- sample
	}
	
	currentdn <- dimnames(newoutput$mcmc[[1]])
	iternames <- seq((burnin+1), (burnin+(sample*thin))-(thin-1), length.out=sample)

	for(i in 1:length(newoutput$mcmc)){
		useits <- iternames
		if(dim(newoutput$mcmc[[i]])[1]!=sample){
			warning(paste('The dimensions of the returned MCMC object did not match the expected sample length for chain ', i, ' - the iteration numbers returned may not be correct', sep=''))
			useits <- useits[1:dim(newoutput$mcmc[[i]])[1]]
		}
		dimnames(newoutput$mcmc[[i]]) <- list(iternames, currentdn[[2]])
	}
	if(class(newoutput$pd)=="mcmc" && !is.na(newoutput$pd)){
		dimnames(newoutput$pd) <- list(iternames, dimnames(newoutput$pd)[[2]])			
	}	

	# Combine runjags objects if necessary:
	if(runjags.object$combine && runjags.object$sample>0){
		burnin <- runjags.object$oldburnin
		
		s <- try(combpd <- if('full.pd'%in%runjags.object$monitor) combine.mcmc(list(runjags.object$pd, newoutput$pd), collapse.chains=FALSE) else NA)
		if(class(s)=='try-error'){
			warning('An unexpected error occured while trying to combine the full.pD from the old and new simulations - it has been removed', call.=FALSE)
			runjags.object$monitor <- runjags.object$monitor[runjags.object$monitor!='full.pd']
			combpd <- NA
		}
		
		combinedoutput <- list(mcmc=combine.mcmc(list(runjags.object$mcmc, newoutput$mcmc), collapse.chains=FALSE), deviance.table=weightedaverage(runjags.object$deviance.table, newoutput$deviance.table, niter(runjags.object$mcmc), niter(newoutput$mcmc)), deviance.sum=weightedaverage(runjags.object$deviance.sum, newoutput$deviance.sum, niter(runjags.object$mcmc), niter(newoutput$mcmc)), pd=combpd, end.state=newoutput$end.state, samplers=newoutput$samplers)

	}else{
		combinedoutput <- list(mcmc=newoutput$mcmc, deviance.table=newoutput$deviance.table, deviance.sum=newoutput$deviance.sum, pd=if('full.pd'%in%runjags.object$monitor) newoutput$pd else NA, end.state=newoutput$end.state, samplers=newoutput$samplers)
		# Never change thin - it messes everything up. Iteration labels will be slightly wrong, that's all.
		# runjags.object$thin <- thin
		# This means if extended, the thin will be what we had here
	}

	# Save some RAM:
	rm(newoutput)
	gcinfo <- gc(FALSE)
	
	# Takes into account (1) the time taken to import, (2) the time taken for the previous sims if any, (3) the simulation run time NOT including any gap between finishing and importing
	timetaken <- (difftime(Sys.time(), starttime) + runjags.object$timetaken + simtimetaken)
	
	expected <- sample
	if(runjags.object$combine) expected <- expected + runjags.object$sample
	if(expected != niter(combinedoutput$mcmc)) warning(paste('Unexpected discrepancy in sample size: expected ', expected, ' iterations but there are actually ', niter(combinedoutput$mcmc), sep=''))
	
	# Call function to calculate summary statistics and plots etc:	
	combinedoutput <- makerunjagsobject(combinedoutput, summarise=runjags.object$summarise, summaryargs=summaryargs, burnin=burnin, sample=sample, thin=thin, model=runjags.object$model, data=runjags.object$data, monitor=runjags.object$monitor, noread.monitor=runjags.object$noread.monitor, modules=runjags.object$modules, factories=runjags.object$factories, response=runjags.object$response, residual=runjags.object$residual, fitted=runjags.object$fitted, method=runjags.object$method, method.options=runjags.object$method.options, timetaken=timetaken)

	stopifnot(class(combinedoutput$end.state)=='runjagsinits')
	
	swcat("Finished running the simulation\n")
	
	# Reset deleting option to what it should be:
	runjags.object$keep.jags.files <- reallydelete
	
#	if(identical(combinedoutput$samplers, NA) && runjags.getOption('nodata.warning'))
#		warning('JAGS is reporting that no samplers were used - ensure that any data has been passed to JAGS correctly')
#  This can happen when just forward sampling from the prior
	
	return(combinedoutput)

}

results.JAGS <- results.jags

