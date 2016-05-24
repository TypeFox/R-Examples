#' @title The runjags class and available S3 methods
#' @name runjags-class
#' @aliases runjags-class runjagsclass runjagsstudy-class runjagsstudyclass as.jags as.runjags as.mcmc.runjags as.mcmc.list.runjags cleanup.jags cleanup.JAGS failed.jags failed.JAGS failedjags fitted.runjags residuals.runjags predict.runjags

#' @export

#' @description
#' Objects of class 'runjags' are produced by \code{\link{run.jags}}, \code{\link{results.jags}} and \code{\link{autorun.jags}}, and contain the MCMC chains as well as all information required to extend the simulation.  These are a number of utility functions associated with these objects.

#' @details 
#' The functions and methods detailed here permit conversion of runjags objects to MCMC objects and to/from jags models created by \code{\link[rjags]{jags.model}}.  There are also S3 methods for print, summary and plot available for runjags class objects - see \code{\link{add.summary}} for details of the arguments available to these.

#' The 'failed.jags' function allows the user to interrogate the details of JAGS models that failed to compile or produce MCMC output.  By default, any simulation folders for models that failed to import are kept until the R session is ended - in some circumstances it may be possible to partially recover the results using \code{\link{results.jags}}.  The cleanup.jags function can be used to remove simulation folders created in the current R session, and is called when the runjags package is unloaded.

#' @examples
#' if(require('rjags')){
#' # Coercion between jags and runjags objects (requires loading the rjags package):
#' data(LINE)
#' jags.model <- LINE
#' runjags.model <- as.runjags(jags.model, monitor=c('alpha','beta'))
#' runjags.model <- extend.jags(runjags.model, method='interruptible')
#' jags.model <- as.jags(runjags.model)
#' # Coercion to MCMC (requires loading the coda package):
#' library('coda')
#' mcmc <- as.mcmc.list(runjags.model)
#' summary(mcmc)
#' }

#' @keywords models

#' @seealso
#' \code{\link{add.summary}} for details on plot, print and summary methods for runjags class objects, \code{\link{extract.runjags}} for a method to extract peripheral information from runjags objects, \code{\link{runjags.options}} for general options available, and \code{\link{run.jags}} and \code{\link{autorun.jags}} for the functions that create objects of this class.

#' @param x an object of class runjags.
#' @param object an object of class runjags.
#' 
#' @param vars an optional character vector of variable names to extract.  If supplied, only variable names in the object supplied with a partial match to anything in 'vars' will be summarised/plotted/extracted.  Note that regular expressions are not allowed, but the caret (^) token can be used to specify the match at the start of a variable name, and a quoted vars will be matched exactly.  Default NA meaning all variables available are returned.
#' 
#' @param add.mutate option to use the inbuild mutate function to produce additional MCMC variables before returning the MCMC object.
#'
#' @param adapt as for \code{\link[rjags]{jags.model}}
#'
#' @param quiet as for \code{\link[rjags]{jags.model}}
#'
#' @param jags.model a model produced by \code{\link[rjags]{jags.model}}

#' @param monitor a character vector of the names of variables to monitor, as for \code{\link{run.jags}}

#' @param modules a character vector of external modules to be loaded into JAGS, either as the module name on its own or as the module name and status separated by a space, for example 'glm on'.  

#' @param factories a character vector of factory modules to be loaded into JAGS.  Factories should be provided in the format '<facname> <factype> <status>' (where status is optional), for example: factories='mix::TemperedMix sampler on'.  You must also ensure that any required modules are also specified (in this case 'mix'). 

#' @param jags the system call or path for activating JAGS.  Default uses the option given in \code{\link{runjags.options}}.

#' @param mutate either a function or a list with first element a function and remaining elements arguments to this function that can be used to add variables to the model output.  See \code{\link{add.summary}} for more details.

#' @param check option to check that the model can be (re)-compiled.

#' @param all.folders option to remove ALL simulation folders created using keep.jags.files=TRUE and not just unsuccessful simulations.

#' @param silent option to suppress feedback when deleting simulation folders.

#' @param show which parts of the failed JAGS simulation to display - options are:  'model', 'data', 'inits', 'output', 'end.state', 'all'

#' @param variable the name of the variable within the JAGS simulation that denotes the residual/fitted variable.  This must be specified to be able to use the residuals and fitted methods.

#' @param show.summary option to show the full summary statistics of the returned models before extracting just the residuals/fitted variable information.

#' @param output the type of output required for the residuals and fitted methods - options are:  'mean', 'mcmc', 'hpd', 'summary', 'runjags'.

#' @param ... additional options to be passed to default methods or additional functions.
NULL
		
#' @rdname runjags-class
#' @method as.mcmc runjags
as.mcmc.runjags <- function(x, vars=NA, add.mutate=TRUE, ...){
	
	m <- as.mcmc.list(x, vars=vars, add.mutate=add.mutate, ...)
	
	# as.mcmc throws an error if there are more than 1 chain - should this be the correct behaviour?  Or stop and say use combine.mcmc???
	if(length(m)>1)
		warning(paste("Combining the ", length(m), " mcmc chains together", sep=""))
	
	m <- as.mcmc(combine.mcmc(m, collapse.chains=TRUE))

	return(m)
	
}


#' @rdname runjags-class
#' @method as.mcmc.list runjags
as.mcmc.list.runjags <- function(x, vars=NA, add.mutate=TRUE, ...){
	
	if(add.mutate) 
		m <- as.mcmc.list(addmutated(x$mcmc, x$summary.pars$mutate))
	else
		m <- as.mcmc.list(x$mcmc)

	selected <- matchvars(checkvalidmonitorname(vars),  varnames(m))
	thevarnames <- dimnames(m[[1]])
	
	m <- as.mcmc.list(lapply(m, function(x) return(x[,selected,drop=FALSE])))
	for(i in 1:length(m)){
		dimnames(m[[i]]) <- list(thevarnames[[1]], thevarnames[[2]][selected])
	}
	
	return(as.mcmc.list(m))
	
}

# # A possible solution to the S3 method issues with as.mcmc.list:
# as.mcmclist.runjags <- function(x, vars=NA, add.mutate=TRUE, ...){
#   
#   if(add.mutate) 
#     m <- as.mcmclist(addmutated(x$mcmc, x$summary.pars$mutate))
#   else
#     m <- as.mcmclist(x$mcmc)
#   
#   selected <- matchvars(checkvalidmonitorname(vars),  varnames(m))
#   thevarnames <- dimnames(m[[1]])
#   m <- m[,varnames(m)[selected],drop=FALSE]	
#   for(i in 1:length(m)){
#     dimnames(m[[i]]) <- list(thevarnames[[1]], thevarnames[[2]][selected])
#   }
#   
#   return(as.mcmclist(m))
#   
# }
# as.mcmclist <- function(x, ...){
#   UseMethod("as.mcmclist")
# }
# as.mcmclist.default <- function(x, ...){
#   return(coda::as.mcmc.list(x, ...))
# }


#' @rdname runjags-class
#' @method as.jags runjags
as.jags.runjags <- function(x, adapt=1000, quiet=FALSE, ...){
		
	passed <- list(...)
  if(length(passed)>0)
    stop(paste('unused argument(s) ', paste(names(passed),collapse=', '), sep=' '))
  
	if(!loadandcheckrjags(FALSE))
		stop('The rjags package is required for jags/runjags conversion tools')
	runjags.object <- x
	
	runjags.object$factories <- checkmodfact(runjags.object$factories, 'factory')
	runjags.object$modules <- checkmodfact(runjags.object$modules, 'module')
	
  if(any(grepl('lecuyer',runjags.object$end.state)) && !any(unlist(runjags.object$modules)=='lecuyer')){
    stop('The lecuyer module is required to compile this model but has not been added to the modules for the runjags object', call.=FALSE)
  }
  
	# Module loading MUST be done before model compilation - doesn't do any harm to re-load modules:
	if(!identical(runjags.object$modules,"")) for(i in 1:length(runjags.object$modules)){
		if(runjags.object$modules[[i]][1]=="runjags"){
			if(runjags.object$modules[[i]][2]=='TRUE'){
				success <- load.runjagsmodule()
			}else{
				success <- unload.runjagsmodule()
			}
		}else{
			if(runjags.object$modules[[i]][2]=='TRUE'){
				success <- try(rjags::load.module(runjags.object$modules[[i]][1]))
				if(is.null(success)) success <- TRUE
			}else{
				suppressWarnings(success <- try(rjags::unload.module(runjags.object$modules[[i]][1])))
				if(is.null(success)) success <- TRUE
			}
		}
		if(inherits(success, 'try-error')) stop(paste("Failed to ", if(runjags.object$modules[[i]][2]=='FALSE') "un", "load the module '", runjags.object$modules[[i]][1], "'", sep=""))			
	}		
	if(!identical(runjags.object$factories,"")) for(i in 1:length(runjags.object$factories)){
		fa <- ""
		try(fa <- as.character(rjags::list.factories(runjags.object$factories[[i]][2])$factory))
		if(!runjags.object$factories[[i]][1] %in% fa) stop(paste("The factory '", runjags.object$factories[[i]][1], "' of type '", runjags.object$factories[[i]][2], "' is not available - ensure any required modules are also provided", sep=""))

		success <- try(rjags::set.factory(runjags.object$factories[[i]][1], runjags.object$factories[[i]][2], as.logical(runjags.object$factories[[i]][3])))
		if(is.null(success)) success <- TRUE

		if(inherits(success, 'try-error')) stop(paste("Failed to ", if(runjags.object$factories[[i]][3]=='FALSE') "un", "set the factory '", runjags.object$factories[[i]][1], "' of type '", runjags.object$factories[[i]][2], "'", sep=""))			
	}


	# See if the model already exists on the runjags object:
	if( 'rjags' %in% names(runjags.object$method.options)){
		
		jags.object <- runjags.object$method.options$rjags
		
	}else{
		
		if(!quiet) swcat("Compiling rjags model", if(adapt>0) paste(" and adapting for ", adapt, " iterations", sep=""), "...\n",sep="")
		flush.console()
	
		model <- textConnection(runjags.object$model)
		dataenv <- list.format(as.character(runjags.object$data))
		inits <- lapply(runjags.object$end.state,list.format)
		o <- capture.output({s <- try({
		if(length(inits[[1]])==0){
			if(as.character(runjags.object$data)==""){
				jags.object <- rjags::jags.model(model, n.chains=length(runjags.object$end.state), n.adapt=adapt, quiet=quiet)
			}else{
				jags.object <- rjags::jags.model(model, data=dataenv, n.chains=length(runjags.object$end.state), n.adapt=adapt, quiet=quiet)				
			}
		}else{
			if(as.character(runjags.object$data)==""){
				jags.object <- rjags::jags.model(model, inits=inits, n.chains=length(runjags.object$end.state), n.adapt=adapt, quiet=quiet)
			}else{
				jags.object <- rjags::jags.model(model, data=dataenv, inits=inits, n.chains=length(runjags.object$end.state), n.adapt=adapt, quiet=quiet)				
			}
		}}, silent=TRUE)})
	
		close(model)
		flush.console()
	
		if(class(s)=="try-error"){
			jagsout <- as.character(s)
			class(jagsout) <- "rjagsoutput"
			assign("model", runjags.object$model, envir=failedjags)
			assign("inits", runjags.object$end.state, envir=failedjags)
			assign("data", runjags.object$data, envir=failedjags)
			assign("output", jagsout, envir=failedjags)
			stop(paste("The following error occured when compiling and adapting the model using rjags:\n ",jagsout,"\nIt may help to use failed.jags(c('model','data','inits')) to see model/data/inits syntax with line numbers",sep=""),call.=FALSE)
		}
		
	}	
	
	# Now check it has compiled:
	
	checkcompiled <- try(stats::coef(jags.object),silent=TRUE)
	
	if(class(checkcompiled)=="try-error"){
		if(!quiet) swcat("Re-compiling rjags model", if(adapt > 0) " and adapting", "...\n", sep="")
		o <- capture.output(s <- try(jags.object$recompile(),silent=TRUE))
		if(class(s)=="try-error"){
			jagsout <- as.character(s)
			class(jagsout) <- "rjagsoutput"
			assign("model", runjags.object$model, envir=failedjags)
			assign("inits", runjags.object$end.state, envir=failedjags)
			assign("data", runjags.object$data, envir=failedjags)
			assign("output", jagsout, envir=failedjags)
			
			if(runjags.getOption('debug')){
				cat('Got the following error - is this to do with RNG name vs type again?\n', jagsout, '\n')
				browser()
			}
			
			stop(paste("The following error occured when re-compiling and adapting the model using rjags:\n ",jagsout,"\nIt may help to use failed.jags(c('model','data','inits')) to see model/data/inits syntax with line numbers",sep=""),call.=FALSE)
		}
		
		checkcompiled <- try(stats::coef(jags.object),silent=TRUE)
		if(class(checkcompiled)=="try-error") stop(paste("There was an unexpected error re-compiling this JAGS model:  ", as.character(checkcompiled), sep=""))
	}
	
	return(jags.object)
}



#' @rdname runjags-class
#' @method as.runjags jags
as.runjags.jags <- function(jags.model, monitor = stop("No monitored variables supplied"), modules=runjags.getOption('modules'), factories=runjags.getOption('factories'), jags = runjags.getOption('jagspath'),  mutate=NA, check = TRUE, ...){
	
	if(!loadandcheckrjags(FALSE))
		stop('The rjags package is required for jags/runjags conversion tools')

	jags.object <- jags.model
	model <- paste(jags.object$model(),collapse="\n")
	data <- dump.format(jags.object$data())
	end.state <- sapply(jags.object$state(internal=TRUE), dump.format)
	
	# To be added when they are added for extend.jags etc:
	response <- NA
	fitted <- NA
	residual <- NA 
	
	# Force eval:
	monitor <- force(monitor)
	
	rjo <- setup.jags(model=model, monitor = monitor, data=data,  n.chains=length(end.state), inits = end.state, modules=modules, factories=factories, response=response, fitted=fitted, residual=residual, jags = jags, method=if(check) "simple" else "rjags", mutate=mutate)
	
	timetaken <- as.difftime(0,units='secs')
	
	blankmcmc <- as.mcmc.list(lapply(1:length(end.state), function(x) return(as.mcmc(matrix(NA,ncol=1,nrow=0)))))			
	
	combinedoutput <- makerunjagsobject(list(mcmc=blankmcmc, deviance.table=NA, deviance.sum=NA, pd=NA, end.state=end.state, samplers=NA, end.time=Sys.time()), summarise=FALSE, summaryargs=list(), burnin=0, sample=0, thin=1, model=rjo$model, data=rjo$data, monitor=monitor, noread.monitor='', modules=rjo$modules, factories=rjo$factories, response=rjo$response, residual=rjo$residual, fitted=rjo$fitted, method=rjo$method, method.options=rjo$method.options, timetaken=timetaken, silent=TRUE)
	
	return(combinedoutput)
	
}

#' @rdname runjags-class
is.runjags <- function(x){
  return(class(x) %in% c("runjags","runjagsbginfo"))
}

#' @rdname runjags-class
cleanup.jags <- function(all.folders=FALSE, silent=FALSE){
	if(all.folders){
		tobin <- c(runjagsprivate$simfolders, runjagsprivate$failedsimfolders)
	}else{
		tobin <- runjagsprivate$failedsimfolders
	}
	
	tobin <- tobin[file.exists(tobin)]
	for(f in tobin){
		unlink(f, recursive=TRUE)
	}

	runjagsprivate$failedsimfolders <- character(0)
	if(length(tobin)==0){
		if(!silent) swcat('No JAGS simulation folders were found\n')		
	}else{
		if(all.folders){
			if(!silent) swcat('Your', length(tobin), 'JAGS simulation folder(s) from this R session were deleted\n')
			runjagsprivate$simfolders <- character(0)
		}else{
			if(!silent) swcat('Your', length(tobin), 'failed JAGS simulation folder(s) from this R session were deleted\n')
		}
	}

	invisible(tobin)
}
cleanup.JAGS <- cleanup.jags


#' @rdname runjags-class
failed.jags <- function(show=c('model','output')){	
	type <- tolower(show)
	possibilities <- c("model", "data", "inits", "output", "end.state")
	if(any(type=='all'))
		type <- possibilities
	if(!all(type %in% possibilities))
		stop(paste('Unsupported failed jags type(s) ', paste(type[!type %in% possibilities],collapse=', '), ' - options are ', paste(possibilities,collapse=', '), sep=''))	
	
	output <- lapply(type, function(x){
		toret <- switch(x, 'model'=failedjags$model, 'data'=failedjags$data, 'inits'=failedjags$inits, 'output'=failedjags$output, 'end.state'=failedjags$end.state)    
	})
	
  	names(output) <- type
	class(output) <- 'failedjags'
		
	return(output)
}
failed.JAGS <- failed.jags


# Not documented specifically:
as.jags <- function(x, adapt=1000, quiet=FALSE, ...){
	UseMethod("as.jags")
}
as.jags.default <- function(x, ...){
	stop("Conversion to jags objects is only possible for specific classes of object (such as objects of class 'runjags')")
}
as.runjags <- function(jags.model, ...){
	UseMethod("as.runjags")
}
as.runjags.default <- function(jags.model, ...){
	stop("Conversion to runjags objects is only possible for specific classes of object (such as objects of class 'jags')")
}




#' @rdname runjags-class
#' @method residuals runjags
residuals.runjags <- function(object, variable=object$residual, show.summary=FALSE, output='mean', ...){
	
	if(length(output)!=1 || class(output)!='character')
		stop('The output specification must be a single character string')
	output <- tolower(output)
	possibilities <- c('mean','mcmc','hpd','summary','runjags')
	if(!output %in% possibilities)
		stop(paste('Unsupported output option "', output, '" - options are ', paste(possibilities,collapse=', '), sep=''))	
	
	object <- checkvalidrunjagsobject(object)
	
	if(is.na(variable)){
		if(is.na(object$response) || is.na(object$fitted))
			stop('No residual variable (or response and fitted variables) specified')
		
		tomonitor <- checkvalidmonitorname(c(object$response, object$fitted))
				
	}else{
		tomonitor <- checkvalidmonitorname(variable)
	}
	
	swcat('Extending the simulation to extract the residuals...\n')
	
	s <- try({
	mcmcout <- extend.jags(object, add.monitor=tomonitor, drop.monitor=object$monitor, combine=FALSE, summarise=FALSE, ...)
	})
	if(class(s)=='try-error')
		stop('An unexpected error occuring while obtaining the MCMC samples of the residuals - consult the error message above to diagnose the problem')
	
	# Check a couple of things if we have calculated the residual:
	if(length(tomonitor)==2){
		  s <- try(respselected <- matchvars(object$response, dimnames(mcmcout$mcmc[[1]])[[2]]), silent=TRUE)
		  if(class(s)=='try-error')
		 	 stop('The response variable specified could not be monitored')
		  s <- try(fittedselected <- matchvars(object$fitted, dimnames(mcmcout$mcmc[[1]])[[2]]), silent=TRUE)
		  	if(class(s)=='try-error')
		  	  stop('The fitted variable specified could not be monitored')
  		  if(length(respselected)!=length(fittedselected))
  		    stop('The length of the response vector did not match the length of the fitted vector')
  		
		nonstochastic <- normalise.mcmcfun(mcmcout$mcmc, normalise = FALSE, warn=FALSE, remove.nonstochastic = FALSE)$nonstochastic
		if(!all(nonstochastic[dimnames(mcmcout$mcmc)[[2]][respselected]]))
			warning('One or more of the response variables is stochastic - the residual calculation may not be correct!')
		
		# A mutate function to post-calculate the residual:
		mutatefun <- list(function(mcmc, response, fitted){
  
		  respselected <- matchvars(response, dimnames(mcmc)[[2]])
		  fittedselected <- matchvars(fitted, dimnames(mcmc)[[2]])
  
		  resid <- mcmc[,respselected,drop=FALSE] - mcmc[,fittedselected,drop=FALSE]
		  newdn <- gsub(response, 'residual', dimnames(mcmc)[[2]][respselected], fixed=TRUE)
		  dimnames(resid) <- list(dimnames(mcmc)[[1]], newdn)
		  return(as.mcmc(resid))
  
		}, response=object$response, fitted=object$fitted)
		
		variable <- 'residual'
	}else{
		mutatefun <- NULL
	}
	
	swcat('Calculating summary statistics...\n')
	s <- try({
		mcmcout <- add.summary(mcmcout, mutate=mutatefun, plots=FALSE, vars=variable, silent.jags=TRUE)
		})
	if(class(s)=='try-error')
		stop('An unexpected error occuring while calculating summary statistics for the residuals - consult the error message above to diagnose the problem')	
	
	stochastic <- !normalise.mcmcfun(as.mcmc.list(mcmcout, add.mutate=TRUE, vars=variable), normalise = FALSE, warn=FALSE, remove.nonstochastic = FALSE)$nonstochastic
	if(!all(stochastic))
		warning('One or more of the residuals is non-stochastic!')
	
	if(show.summary)
		print(mcmcout)
	
	summarystats <- summary(mcmcout)
	lu <- grepl('Lower',dimnames(summarystats)[[2]]) | grepl('Upper',dimnames(summarystats)[[2]])
	
	swcat('Finished extracting the residuals\n')
	
	if(output=='mean')
		return(summarystats[,'Mean'])
	if(output=='summary')
		return(summarystats)
	if(output=='hpd')
		return(t(summarystats[,lu]))
	
	if(output=='mcmc')
		return(combine.mcmc(mcmcout, collapse.chains=FALSE, add.mutate=TRUE, vars=variable))
	
	if(output=='runjags')
		return(mcmcout)
	
}

#' @rdname runjags-class
#' @method fitted runjags

fitted.runjags <- function(object, variable=object$fitted, show.summary=FALSE, output='mean', ...){
	
	if(length(output)!=1 || class(output)!='character')
		stop('The output specification must be a single character string')
	output <- tolower(output)
	possibilities <- c('mean','mcmc','hpd','summary','runjags')
	if(!output %in% possibilities)
		stop(paste('Unsupported output option "', output, '" - options are ', paste(possibilities,collapse=', '), sep=''))	
	
	object <- checkvalidrunjagsobject(object)
	
	if(is.na(variable))
		stop('No response variable specified')
		
	tomonitor <- checkvalidmonitorname(variable)
	
	swcat('Extending the simulation to extract the fitted variable monitor...\n')
	
	s <- try({
	mcmcout <- extend.jags(object, add.monitor=tomonitor, drop.monitor=object$monitor, combine=FALSE, summarise=FALSE, ...)
	})
	if(class(s)=='try-error')
		stop('An unexpected error occuring while obtaining the MCMC samples of the fitted variable - consult the error message above to diagnose the problem')
	
	swcat('Calculating summary statistics...\n')
	
	s <- try({
		mcmcout <- add.summary(mcmcout, mutate=NULL, plots=FALSE, vars=variable, silent.jags=TRUE)
		})
	if(class(s)=='try-error')
		stop('An unexpected error occuring while calculating summary statistics for the fitted variable - consult the error message above to diagnose the problem')	
	
	stochastic <- !normalise.mcmcfun(as.mcmc.list(mcmcout, add.mutate=TRUE, vars=variable), normalise = FALSE, warn=FALSE, remove.nonstochastic = FALSE)$nonstochastic
	if(!all(stochastic))
		warning('One or more of the fitted variables is non-stochastic!')
	
	if(show.summary)
		print(mcmcout)
	
	summarystats <- summary(mcmcout)
	lu <- grepl('Lower',dimnames(summarystats)[[2]]) | grepl('Upper',dimnames(summarystats)[[2]])
	
	swcat('Finished extracting the fitted variable\n')
	
	if(output=='mean')
		return(summarystats[,'Mean'])
	if(output=='summary')
		return(summarystats)
	if(output=='hpd')
		return(t(summarystats[,lu]))
	
	if(output=='mcmc')
		return(combine.mcmc(mcmcout, collapse.chains=FALSE, add.mutate=TRUE, vars=variable))
	
	if(output=='runjags')
		return(mcmcout)
	
}

predict.runjags <- function(object, newdata=NULL, ...){
	
	stop('Not implemented yet')
	
	# If N is in the fitted model data but not in newdata, set N=nrow(newdata)
	# Remove existing monitors, set response to monitor
	# Take 1 random row or mean of estimates as data (variable)
	# Return options as for resid/fitted but plus 'random'
	# Can have type = c("link", "response") (monitor fitted or response) - need to know response to remove either way
	# Set sampler.warning = FALSE and warn if the samplers are used
	# Chat in help file about what is and isn't integrated over wrt coefficients	
	# Number of grouping variables may have to be changed as well - can be too many in the original model without problems, so include all possible levels of factors then. 
	# All of this will require specified coefs to be removed from inits, and an adaptation and burnin might be required if there are still likelihoods being used - I can check if samplers is NA and give a warning though
	
}