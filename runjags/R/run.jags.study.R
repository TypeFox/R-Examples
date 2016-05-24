#' @title Drop-k and simulated dataset studies using JAGS
#' @name run.jags.study
#' @aliases run.jags.study run.JAGS.study drop.k drop.k.jags drop.k.JAGS
#' @export

#' @description
#' These functions can be used to fit a user specified JAGS model to multiple datasets with automatic control of run length and convergence, over a distributed computing cluster such as that provided by snow.  The results for monitored variables are compared to the target values provided and a summary of the model performance is returned.  This may be used to facilitate model validation using simulated data, or to assess model fit using a 'drop-k' type cross validation study where one or more data points are removed in turn and the model's ability to predict that datapoint is assessed. 
#'

#' @details
#' The drop.k function is a wrapper to run.jags.study for the common application of drop-k cross validation studies on fitted JAGS models.  The run.jags.study function is more flexible, and can be used for validating the performance of a model against simulated data with known parameters.  For the latter, a user-specified function to generate suitable datasets to analyse is required.

#' @keywords methods

#' @return
#' An object of class \code{\link{runjagsstudy-class}}, containing a summary of the performance of the model with regards to the target variables specified.  If record.chains=TRUE, an element named 'runjags' containing a list of all the runjags objects returned will also be present.  Any error messages given by individual simulations will be contained in the $errors element of the returned list.

#' @references
#' M. J. Denwood, "runjags: An R Package Providing Interface Utilities, Distributed Computing Methods and Additional Distributions For MCMC Models in JAGS," Journal of Statistical Software, [Under review].

#' @seealso
#' \code{\link{autorun.jags}} for the underlying methods used to run simulations to convergence, and \code{\link{runjagsstudy-class}} for details of the returned object

#' @examples
#' # For examples of usage see the following vignette:
#' \dontrun{
#' vignette('userguide', package='runjags')
#' }

#' @param runjags.object an object of class \code{\link{runjagsstudy-class}} on which to perform the drop-k analysis

#' @param dropvars the variable(s) to be eliminated from the data so that the ability of the model to predict these datapoints can be assessed.  The variable can be specified as a vector, or as a single character for which partial matching will be done.  Array indices can be used, but must be specified as a complete range e.g. variable[2:5,2] is permitted, but variable[,2] is not because the first index is empty

#' @param k the number of datapoints to be dropped from each individual simulation.  The default of 1 is a drop-1 study (also called a leave-one-out cross validation study).

#' @param simulations the number of datasets to run the model on.  For drop.k the default is to use the number of unique datapoints, resulting in a drop-1 study.  If the specified number of simulations is different to the number of unique datapoints, the datapoints are dropped randomly between simulations.

#' @param model the JAGS model to use, in the same format as would be specified to \code{\link{run.jags}}.

#' @param datafunction a function that will be used to specify the data. This must take either zero arguments, or one argument representing the simulation number, and return either a named list or character vector in the R dump format containing the data specific to that simulation.  It is possible to specify any data that does not change for each simulation using a #data# <variable> tag in  the model code.

#' @param targets a named list of variables (which can include vectors/arrays) with values to which the model outputs are compared (if stochastic).  The target variable names are also automatically included as monitored variables.  

#' @param confidence a probability (or vector of probabilities) to use when calculating the proportion of credible intervals containing the true target value.  Default 95\% CI.  

#' @param record.chains option to return the full runjags objects returned from each simulation as a list item named 'runjags'.  

#' @param max.time the maximum time for which each individual simulation is allowed to run by the underling autorun.jags function. Acceptable units include 'seconds', 'minutes', 'hours', 'days', 'weeks', or the first letter(s) of each.  Default is 15 minutes.  

#' @param silent.jags option to suppress all JAGS output, even for simulations run locally.  If set to FALSE, there is no guarantee that the output will be displayed in sequential order between the parallel simulations.  Default TRUE.

#' @param n.cores the maximum number of cores to use for parallel simulations.  Default value uses \code{\link[parallel]{detectCores}}, or a minumum of 2.  Ignored if cl is supplied, or if parallel.method does not take a cl argument.

#' @param parallel.method a function that will be used to call the repeated simulations.  This must take the first two arguments 'X' and 'FUN' as for \code{\link[base]{lapply}}, with other optional arguments passed through from the parent function call.  Default uses \code{\link[parallel]{parLapply}}, but \code{\link[base]{lapply}} or \code{\link[parallel]{mclapply}} could also be used.  

#' @param export.cluster a character vector naming objects to be retrieved from the parent frame of the function call and made available to the cluster nodes.  This may be useful if the initial values specified for the model are required to be extracted from the working environment, however it may be preferable to specify a function for inits instead.

#' @param inits as for \code{\link{run.jags}}, except that it is not permitted to be an environment.  It is recommended to a function to return appropriate initial values (which may depend on the data visible when the function is evaluated).

#' @param ... optional arguments to be passed to \code{\link{autorun.jags}}, or to the parallel method function (such as 'cl'). 
NULL


#' @rdname run.jags.study
drop.k <- function(runjags.object, dropvars, k=1, simulations=NA, ...){
	
	runjags.object <- checkvalidrunjagsobject(runjags.object)
	if(runjags.object$sample==0)
		stop('Cannot perform a drop-k study on a simulation with no samples')
	
	if(missing(dropvars))
		stop('The name of the variable to be dropped must be specified as dropvars')
	
	# Some options that aren't allowed:
	if(any(c('model','datafunction','inits', 'monitor', 'data', 'n.chains', 'combine', 'datalist', 'initlist', 'export.cluster', 'test') %in% names(list(...)))){
		stop("The following options are not permitted for drop.k: 'model','datafunction','inits', 'monitor', 'data', 'n.chains', 'combine', 'datalist', 'initlist', 'export.cluster', 'test'")
	}
	
	# Harvest required info:
	model <- runjags.object$model
	n.chains <- length(runjags.object$end.state)
	
	# These will be used quite a lot - arrays are explicitly rolled out to make overwriting easier:
	basedata <- getjagsnames(list.format(runjags.object$data))
	testinits <- lapply(runjags.object$end.state, function(x) return(getjagsnames(list.format(x))))
	baseinits <- lapply(runjags.object$end.state, list.format)
	getarraynames <- getarraynames
	getjagsnames <- getjagsnames
	
	# Match the data to vars:
	jagsnames <- names(basedata)
	vars <- jagsnames[matchvars(checkvalidmonitorname(dropvars), jagsnames)]
	newinits <- basedata[vars]
	
	if(k > length(vars)){
		stop('Specified "k" was greater than the number of data points', call.=FALSE)
	}
	if(k < 1)
		stop('Specified "k" was less than 1', call.=FALSE)

	# One way to specify drop all:
	if(k == length(vars)){
		if(!is.na(simulations))
			warning('Value for "simulations" ignored as k == length(vars)', call.=FALSE)	
		simulations <- 1
	}
	
	if(length(simulations)!=1 || is.na(simulations)){
		if(k!=1)
			stop('A single positive integer must be provided for "simulations" if k>=2', call.=FALSE)		
		simulations <- length(vars)
	}
	
	# Other way to specify drop all:
	if(simulations==1){
		if(k!=1 && k!=length(vars)){
			warning('Value for "k" ignored as "simulations" was specified as 1', call.=FALSE)
		}
		k <- length(vars)
	}
	
	# The inits function takes base inits and base data and data removed and acts appropriately
	
	# datafunction should be able to see basedata and baseinits and vars as this is the parent frame
	dfdone <- FALSE
	avs <- FALSE
	if(k==1 && simulations==length(vars)){
		# Non random drop-1:
		datafunction <- function(i){
			thisdata <- basedata
			thisdata[vars[i]] <- as.numeric(NA)
			return(getarraynames(thisdata))
		}
		dfdone <- TRUE
		avs <- FALSE
	}
	
  if(k==length(vars)){
		# Non random drop-all - this is probably not a useful thing to want to do....
		datafunction <- function(i){
			thisdata <- basedata
			thisdata[] <- as.numeric(NA)
			return(getarraynames(thisdata))
		}
		dfdone <- TRUE
		avs <- FALSE
		if(simulations!=1 && !is.na(simulations))
			warning('Ignoring provided value for simulations, as "k" is equal to the number of variables')
		simulations <- 1		
	}
	if(!dfdone){
		# Otherwise, random drop-k:
		datafunction <- function(i){
			thisdata <- basedata
			makena <- sample(1:length(vars), k)
			thisdata[vars[makena]] <- as.numeric(NA)
			return(getarraynames(thisdata))
		}		
		avs <- TRUE
	}
	
    inits <- function(chain){
		# We have the data returned by datafunction visible:
		thisdata <- getjagsnames(data)
		# And the baseinits exported onto the cluster:
		baseinits <- baseinits
		# And the new inits:
		# Taking them from newinits only works if the array is complete - safer to take everything in data and make it NA:
		# thisinits <- newinits
		# makena <- which(!is.na(thisdata[names(thisinits)]))
		# thisinits[makena] <- as.numeric(NA)
		thisinits <- thisdata
		thisinits[] <- as.numeric(NA)
		# Find which are NA in the data:
		bringback <- names(which(is.na(thisdata)))
		thisinits[bringback] <- newinits[bringback]
		# And return all inits:
		toreturn <- c(baseinits[[chain]], getarraynames(thisinits))
		toreturn <- toreturn[sapply(toreturn, function(x) return(!all(is.na(x))))]
		# Needs to be dump format to suppress finding #inits# from the model on the cluster:
		return(dump.format(toreturn))
    }

	results <- run.jags.study(simulations, model=model, datafunction=datafunction, data='', targets=basedata[vars], n.chains=n.chains, inits=inits, export.cluster=list(baseinits=baseinits, newinits=newinits, getjagsnames=getjagsnames, getarraynames=getarraynames), forceaverage=avs, ...)
	
	return(results)
	
}
drop.k.jags <- drop.k
drop.k.JAGS <- drop.k.jags

#' @rdname run.jags.study
run.jags.study <- function(simulations, model, datafunction, targets=list(), confidence=0.95, record.chains=FALSE, max.time="15m", silent.jags=TRUE, parallel.method=parLapply, n.cores=NA, export.cluster=character(0), inits=list(), ...){
	
	# ... is passed either to autorun.jags (add.monitor not allowed, combine not allowed, maybe others) or to parallel.method function
	
	# Reset failedjags stuff:
	failedjags$model <- NA
	failedjags$data <- NA
	failedjags$inits <- NA
	failedjags$output <- NA
	failedjags$end.state <- NA
	
	# Force evaluation:
	parallel.method <- parallel.method	
	simulations <- simulations
	model <- model
	datafunction <- datafunction
	
	passthrough <- list(...)	
	if(any(names(passthrough)=='runjags.options')){
		warning('The runjags.options argument is deprecated - use autorun.jags arguments directly in the run.jags.study function call')
		passthrough <- c(passthrough, passthrough$runjags.options)
		passthrough <- passthrough[names(passthrough)!='runjags.options']
	}
	ncores <- n.cores
	if(identical(ncores, NA)){
		ncores <- max(2, suppressWarnings(parallel::detectCores(TRUE)))  # never less than 2 cores!
	}
	# fromdropk is ust a hack to suppress warnings, forceaverage is set by dropk:
	fromdropk <- FALSE
	forceaverage <- FALSE
	if(any(names(passthrough)=='forceaverage')){
		forceaverage <- passthrough$forceaverage
		passthrough <- passthrough[names(passthrough)!='forceaverage']
		fromdropk <- TRUE
	}	
	# We always need to test so we can check and see which monitors are in the model and which in mutate, and to make sure all targets are found
	# Need to add code to check targets are found etc - at the moment just suppressing warnings
	if(any(names(passthrough)=='test')){
		warning('The test argument is deprecated and was ignored')
	}	
	test <- !fromdropk
	
	# Some checks:
	if(!is.numeric(simulations) || length(simulations)!=1 || simulations!=as.integer(simulations) || simulations<1)
		stop("The value for simulations must be a single positive integer")
	if(!length(confidence)>0 || any(confidence <= 0) || any(confidence > 1))
		stop("The confidence variable must be a numeric vector between 0 and 1")
	
	
	# Set the data scoping - will only be used locally:
	if(!any(names(passthrough)=='data'))
		passthrough$data <- parent.frame()
	
	# The inits should not be an environment as the environment may well not be available on the nodes - can't just make it a list as that would copy everything
	if(!identical(inits, list()) && (class(inits)=='environment' || (class(inits)=='list' && all(sapply(inits,class)=='environment'))))
		stop('The inits argument is not allowed to be an environment for JAGS studies - please use a list or function for inits instead')
	
	# Expected to be in this list:
	passthrough$inits <- inits
	
    if(!fromdropk && !identical(export.cluster, character(0))){
		if(!is.character(export.cluster))
			stop("The argument supplied to 'export.cluster' must be a character vector of variables visible in the working environment")			
 		
		cln <- export.cluster                                                                                                      
        s <- try(objects(name=export.cluster, all.names=TRUE, envir=parent.frame()), silent=FALSE)
		if(class(s)=='try-error')
			stop('One or more variables specified by export.cluster were not found', call.=FALSE)
		export.cluster <- lapply(export.cluster, get, envir=parent.frame())
        names(export.cluster) <- cln                                                                                               
	}else{
		export.cluster <- list()
	}   
	
	# Get the model stuff
	modelsetup <- read.jagsfile(model)
	autodata <- modelsetup$autodata
	maindata <- modelsetup$data	
	if(modelsetup$model=="model{\n\n}\n")
		stop("No valid model was specified or found in the model block")
	##  Get the data
	if(class(passthrough$data)%in%c('runjagsdata','character')){
		if(!fromdropk && (!identical(maindata, NA) || !identical(autodata, NA)) && runjags.getOption('blockignore.warning'))
			warning('Data specified in the model file or using #data# are ignored when a character string is given as the argument to data', call.=FALSE)
		maindata <- NA
		autodata <- NA
		modeldata <- list.format(passthrough$data)
	}else{
		# Unless data is already a character, get the data - ignore inits as they are evaluated on the cluster nodes (except if they are a character, like from drop.k):
		modeldata <- list.format(checkdataformat(arg=passthrough$data, block=maindata, auto=autodata, n.chains=NA, data.type=TRUE, evalscope=NULL)$combined)
	}
	passthrough$data <- NULL
	# We will actually use the fixed data as a character string, which would mean that any #data# lines and blocks in the model will be ignored with a warning ... so I have to remember to temporarily deactivate the warning option later
	
	
	# First check the passed arguments match to something in autorun.jags / parallel.method:
	runjags.args <- getargs(c('autorun.jags','add.summary','parallel.method'), passthrough, returnall=FALSE, 'run.jags.study')
	runjags.args$max.time <- max.time
	
	# Then separate the parallel method options:
	parallel.options <- runjags.args[names(runjags.args) %in% names(formals(parallel.method))]
	runjags.args <- runjags.args[! names(runjags.args) %in% names(parallel.options)]

	# Now remove any passed data list or initlist etc:
	if(any(names(runjags.args)%in%c("datalist","initlist")))
		warning("A datalist and/or initlist supplied was ignored - these arguments are deprecated")
	runjags.args <- runjags.args[!names(runjags.args)%in%c("datalist","initlist")]
	
	if(!silent.jags && !identical(parallel.method, lapply)){
		if(!runjags.getOption('silent.runjags')){
			warning("The progress of JAGS simulations on a remote cluster may not be correctly displayed!", call.=FALSE)
		}
#		silent.jags <- TRUE
	}
	runjags.args$silent.jags <- silent.jags
	runjags.args$batch.jags <- TRUE
	
	if(!any(names(runjags.args)=="method")){
		runjags.args$method <- expression(if(loadandcheckrjags(FALSE, silent=TRUE)) 'rjags' else 'interruptible')
	}else{
		runjags.args$method <- getrunjagsmethod(runjags.args$method)
		if(runjags.args$method%in%runjagsprivate$parallelmethod)
			stop("Cannot use parallel chains for the same simulation with the run.jags.study function")
	}
	# I was having problems with R GUI but I think it was to do with Fork clusters - PSOCK ones work fine
#	if(.Platform$GUI=='AQUA' && runjags.args$method!='rjags')
#		stop('Only the rjags method is supported for JAGS studies using the R.app for Mac OS X - for other methods try using the terminal instead', call.=FALSE)
	
	if(any(names(runjags.args)=='interactive'))
		stop('Cannot use interactive mode for autorun.jags')
	
	if(!any(names(runjags.args)=="summarise")){
		runjags.args$summarise <- FALSE
	}		
	if(any(names(runjags.args)=="plots")){
		if(runjags.args$plots && !runjags.args$summarise){
			warning('Cannot produce plots if summarise=FALSE')
			runjags.args$plots <- FALSE
		}
	}		
	
	runjags.args$model <- model
	
		
#	if(grepl("#inits#",model) && !all(sapply(runjags.args$inits,class)%in%c('character','runjagsinits')) && !fromdropk)
#		warning('The #inits# statement in the model is may give problems on remote machines - to pass initial values to simulation studies safely, use the inits argument (this could be a function, which is evaluated separately for each simulation with "data" in the working environment)', call.=FALSE)
	

	st <- Sys.time()
	swcat(paste("\nStarting a JAGS study at", format(st, format="%H:%M"), "\n"))
	
	# Now check that the datafunction gives us some data and if it does include it with the setup call
	if(!is.null(datafunction)){
		if(!is.function(datafunction) || length(formals(datafunction))>1){
			stop("The datafunction argument provided must be a function with either zero arguments or a single argument representing the simulation number")
		}		

		data <- vector('list',length=simulations)
		alreadywarned <- FALSE
		for(i in 1:simulations){
			if(length(formals(datafunction))==0) thedata <- datafunction() else thedata <- datafunction(i)
			if(class(thedata)=="character") thedata <- list.format(thedata)
			if(class(thedata)!="list") stop("The data function must return either a named list or a character string representing the data for that iteration")
			
			tdata <- modeldata			
			# Check all names are unique and if not over-write tdata (specified in the model) with thedata (specified in the function) with a warning:
			if(any(names(tdata) %in% names(thedata))){
				
				duplicated <- names(tdata)[which(names(tdata) %in% names(thedata))]				
				if(!alreadywarned) warning(paste("The following data variable", if(length(duplicated)>1) "s were" else " was", " specified both in the model and the data function: ", paste(duplicated,collapse=", "), " (the data function values were used)", sep=""))
				alreadywarned <- TRUE
				
				for(d in duplicated) tdata[d] <- thedata[d]
					
				notduplicated <- names(thedata)[!names(thedata) %in% duplicated]
				tdata <- c(tdata, thedata[notduplicated])
				
			}else{
				tdata <- c(tdata, thedata)
			}

			data[[i]] <- dump.format(tdata)
			class(data[[i]]) <- "runjagsdata"
			valid <- checkvalidforjags(data[[i]])	
			if(!valid$valid) stop(paste("The following problem was identified in the data for simulation ", i, ":  ", valid$probstring, sep=""))				
			
		}
		dc <- sample(1:simulations, 1)

		runjags.args$data <- as.character(data[[dc]])
		if(test) swcat("Testing the model and data for simulation ", dc, "...\n",sep="") else swcat("Checking the supplied model definition (and data for simulation ", dc, ")...\n",sep="")		
	}else{
		data <- NULL
		if(test) swcat("Testing the model...\n",sep="") else swcat("Checking the supplied model definition...\n",sep="")		
	}
	
	# Turn list targets into a vector so they can be added to monitored variables:
	if(class(targets)=='list'){
		if(identical(targets,list()))
			stop("A named list or numeric vector of variables and their values must be provided to the 'target' option")
		targets <- getjagsnames(targets)
	}
	if(!class(targets)%in%c("numeric","integer") || length(targets)==0 || any(names(targets)==""))
		stop("The targets variable must be a named list or numeric vector of variable(s) on which to assess the model's performance")
	
	# Now get the full argument list for autorun.jags with these modifications above:
	runjags.args <- getargs(c('autorun.jags','add.summary'), runjags.args, returnall=TRUE)
	donteval <- c('inits','data', 'method', 'jags', 'tempdir')
	for(i in which(!names(runjags.args)%in%donteval)){
		if(!is.null(eval(runjags.args[[i]]))) runjags.args[[i]] <- eval(runjags.args[[i]])
	}

	# If no monitors have been specified, leave it as a zero length character to be added to
	runjags.args$model <- modelsetup$model
	runjags.args$monitor <- c(modelsetup$monitor, eval(runjags.args$monitor))
	if(length(runjags.args$monitor)==1 && is.na(runjags.args$monitor)) runjags.args$monitor <- character(0)
	runjags.args$monitor <- unique(c(runjags.args$monitor, names(targets)))

	# Use a dummy run.jags function call with sample=0 - get rid of add.summary arguments:
	args <- runjags.args[which(names(runjags.args) %in% names(formals(run.jags)))]
	args$summarise <- FALSE
	args$adapt <- 0
	args$burnin <- 0
	args$sample <- 0
	args$silent.jags <- TRUE	
	for(i in which(names(args) %in% c('data', 'method', 'jags', 'tempdir'))){
		if(!is.null(eval(args[[i]]))) args[[i]] <- eval(args[[i]])
	}
	oldmb <- runjags.getOption('blockignore.warning')
	on.exit(runjags.options(blockignore.warning=oldmb))
	
	
	runjags.options(blockignore.warning=FALSE)
	obj <- do.call("run.jags", args=args)
	runjags.options(blockignore.warning=oldmb)
	on.exit()
	
#	if(any(obj$modules=="runjags") && (!eval(runjags.args$method)%in%runjagsprivate$rjagsmethod)) stop("The builtin runjags module is only available for rjags and rjparallel methods")

#  Hope this isn't still required....	
#	if(!runjags.args$summarise && any(obj$monitor=="dic")){
#		runjags.args$summarise <- TRUE
#	}
	
	if(test){
		testoptions <- runjags.args[which(names(runjags.args) %in% names(formals(extend.jags)))]
		testoptions <- c(testoptions, list(runjags.object=obj))
		testoptions$adapt <- 1000
		testoptions$sample <- 1000
		testoptions$combine <- FALSE
		suppressWarnings(testr <- try(do.call("extend.jags", args=testoptions)))   # suppress warnings about invalid monitors
		if(class(testr)=="try-error") stop("The test model returned an error - the simulation study was aborted")
		swcat("The model runs OK\n")
	}else{
		swcat("The model compiles OK\n")
	}
	
	flush.console()
		
	if(is.null(data)){
		DATAS <- lapply(1:simulations, function(x) return(obj$data))
	}else{
		DATAS <- data	
	}
	X <- 1:simulations
	FUN <- function(x, DATAS=DATAS, runjags.args=runjags.args){
		
		runjags.args$method <- eval(runjags.args$method)
		runjags.args$jags <- eval(runjags.args$jags)
		runjags.args$tempdir <- eval(runjags.args$tempdir)
		
		# Detect common problems:
		if(!require("runjags") || package_version(utils::packageDescription('runjags', fields='Version'))<1)
			stop(paste("The runjags package (version >=1.0.0) is not installed (or failed to load) on the cluster node '", Sys.info()['nodename'], "'", sep="")) 
		if(runjags.args$method%in%c('rjags','rjparallel') && !requireNamespace("rjags"))
			stop(paste("The rjags package is not installed (or failed to load) on the cluster node '", Sys.info()['nodename'], "' - try specifying method='simple'", sep="")) 

		if(runjags.args$method%in%c('rjags','rjparallel')){
			
			# Module loading MUST be done before model compilation - doesn't do any harm to re-load modules:
			if(!identical(runjags.args$modules,"")) for(i in 1:length(runjags.args$modules)){
				if(runjags.args$modules[[i]][1]=="runjags"){
					if(runjags.args$modules[[i]][2]=='TRUE'){
						success <- load.runjagsmodule()
					}else{
						success <- unload.runjagsmodule()
					}
				}else{
					if(runjags.args$modules[[i]][2]=='TRUE'){
						success <- try(rjags::load.module(runjags.args$modules[[i]][1]))
						if(is.null(success)) success <- TRUE
					}else{
						suppressWarnings(success <- try(rjags::unload.module(runjags.args$modules[[i]][1])))
						if(is.null(success)) success <- TRUE
					}
				}
				if(inherits(success, 'try-error')) stop(paste("Failed to ", if(runjags.args$modules[[i]][2]=='FALSE') "un", "load the module '", runjags.args$modules[[i]][1], "'", sep=""))
			}		
			if(!identical(runjags.args$factories,"")) for(i in 1:length(runjags.args$factories)){
				fa <- ""
				try(fa <- as.character(rjags::list.factories(runjags.args$factories[[i]][2])$factory))
				if(!runjags.args$factories[[i]][1] %in% fa) stop(paste("The factory '", runjags.args$factories[[i]][1], "' of type '", runjags.args$factories[[i]][2], "' is not available - ensure any required modules are also provided", sep=""))

				success <- try(rjags::set.factory(runjags.args$factories[[i]][1], runjags.args$factories[[i]][2], as.logical(runjags.args$factories[[i]][3])))
				if(is.null(success)) success <- TRUE

				if(inherits(success, 'try-error')) stop(paste("Failed to ", if(runjags.args$factories[[i]][3]=='FALSE') "un", "set the factory '", runjags.args$factories[[i]][1], "' of type '", runjags.args$factories[[i]][2], "'", sep=""))			
			}

		}else{
			if(any(runjags.args$modules=="runjags")) stop("The runjags module is only available using the rjags method; to use the functions provided with other methods install (and specify using the module argument) the 'paretoprior' standalone module")		
		}
		
		thedata <- DATAS[[x]]
		runjags.args$data <- thedata
		
		oldmb <- runjags.getOption('blockignore.warning')
		on.exit(runjags.options(blockignore.warning=oldmb))
		runjags.options(blockignore.warning=FALSE)
		
		if(runjags.args$silent.jags){
			output <- capture.output({
				suppressWarnings(result <- try(do.call("autorun.jags", args=runjags.args)))   # suppress warnings about invalid monitors
				})
		}else{
			suppressWarnings(result <- try(do.call("autorun.jags", args=runjags.args)))    # suppress warnings about invalid monitors
		}
		
		if(!runjags.args$silent.jags)
			swcat("Finished running simulation ", x, " of ", length(DATAS), "\n", sep="")
		
		flush.console()
		
		return(result)
		
	}
	
	if(simulations>1 && all(sapply(2:simulations,function(x) return(identical(DATAS[[1]],DATAS[[x]]))))) warning("The data provided is the same for each simulation")
	
	bigfun <- any(names(formals(parallel.method))=="FUN")
	if(!bigfun && !any(names(formals(parallel.method))=="fun")) stop("The provided parallel apply function does not take the argument 'fun' or 'FUN'")
	if(!any(names(formals(parallel.method))=="X")) stop("The provided parallel apply function does not take the argument 'X'")
	if(any(names(formals(parallel.method)) %in% c("DATAS", "runjags.args"))) stop("The arguments 'DATAS' and 'runjags.args' are used internally and can't be passed to the provided parallel apply function")
	
	swcat("Calling autorun.jags for ", simulations, " simulation", if(simulations>1) "s", "...\n",sep="")
	
	# Now set up the cluster and function call:
	if(identical(parallel.method, parLapply) && !any(names(parallel.options)=="cl")){
		if(is.na(ncores)){
			ncores <- 2
			warning("Unable to detect the available number of cores on your machine - using 2 cores as a default")
		}

		ncores <- min(ncores, simulations)
		parallel.options$cl <- parallel::makeCluster(ncores)    # having problems with Fork cluster on mac gui for some reason - this will always do PSOCK
		on.exit(stopCluster(parallel.options$cl))
	}
	if(any(names(parallel.options)=="cl")){
		clname <- paste(summary(parallel.options$cl)[1,'Class'], capture.output(print(parallel.options$cl)))
		clname <- gsub('forknode socket', 'Fork', clname)
		clname <- gsub('SOCKnode socket', 'PSOCK', clname)			
		swcat('Using a ', clname, '\n', sep='')
		flush.console()
	}
	parallel.options$X <- 1:simulations
	if(bigfun) parallel.options$FUN <- FUN else parallel.options$fun <- FUN 
	parallel.options$DATAS <- DATAS
	parallel.options$runjags.args <- runjags.args
	
	# Export variables - no longer using as inits does the job:
	if(!identical(export.cluster, list()) && !identical(parallel.method, lapply)){
		clusterExport(parallel.options$cl, varlist=names(export.cluster), envir=as.environment(export.cluster))
	}
	
	# And run the function on the cluster:
	success <- try({
		if(runjags.args$silent.jags){
			output <- capture.output({
				results <- do.call('parallel.method', args=parallel.options)
			})
		}else{
			results <- do.call('parallel.method', args=parallel.options)
		}
	})
	flush.console()
	
	if(inherits(success, 'try-error')){
		stop("An unexpected error occured - ensure that the model runs using run.jags on the first dataset, and try using lapply as the parallel.method to debug")
	}
	
	crashed <- sapply(results, class)=="try-error"
	error <- sapply(results[crashed],as.character)

	if(length(error)!=0){
		names(error) <- paste("simulation.",which(crashed),sep="")
		if(eval(runjags.args$method) %in% runjagsprivate$rjagsmethod){
			class(error) <- "rjagsoutput"
		}else{
			class(error) <- "runjagsoutput"
		}
	}
	
	if(all(crashed)){
		failedjags$output <- results
		failedjags$model <- model
		if(any(sapply(results,function(x) return(grepl('not found', x))))){
			if(test) stop("All simulations returned an error (see ?failed.jags) - ensure that the model runs using run.jags on the first dataset, try using lapply as the parallel.method to debug, and remember to use a named list to inits if required for initial values only visible locally") else stop("All simulations returned an error (see ?failed.jags) - ensure that the model runs using run.jags on the first dataset, and remember to use a named list to inits if required for initial values only visible locally")
		}else{
			if(test) stop("All simulations returned an error (see ?failed.jags) - ensure that the model runs using run.jags on the first dataset, and try using lapply as the parallel.method to debug") else stop("All simulations returned an error (see ?failed.jags) - ensure that the model runs using run.jags on the first dataset")
		}
	}

	swcat("Finished running the simulations\n")
	flush.console()

	if(simulations != length(results))
		stop("A different number of results was returned to that specified by simulations - ensure that the model runs using run.jags on the first dataset, and try using lapply as the parallel.method to debug")
	
	# Farm this out into another function so we can easily have a retrieve.jags.study function one day - targets will have to be added to the runjagsstudy.suspended object:
	success <- try({
	retval <- summarise.jags.study(results=results[!crashed], targets=targets, confidence=confidence, separatesingles=!forceaverage)
	})
	if(inherits(success, 'try-error')){
		warning('An unknown error occured while calculating summary statistics')
		retval <- list(simulations=simulations, model=obj$model, targets=targets, monitor=unique(c(obj$monitor,names(targets))), dropk=fromdropk, datafunction=datafunction, data=DATAS, crashed=crashed, errors=error)
	}else{
		retval <- c(retval, list(simulations=simulations, model=obj$model, targets=targets, monitor=unique(c(obj$monitor,names(targets))), dropk=fromdropk, datafunction=datafunction, data=DATAS, crashed=crashed, errors=error))
	}
	if(record.chains) retval <- c(retval, list(runjags=results))
	
	class(retval) <- "runjagsstudy"
	
	et <- Sys.time()
	swcat("Finished summarising results\n")
	
	swcat(paste("Finished JAGS study at ", format(et, format="%H:%M"), " (total time taken:  ", timestring(st, et), ")\n\n", sep=""))
	
	return(retval)
	
}


summarise.jags.study <- function(results, targets, confidence, separatesingles){
	
	if(!length(confidence)>0 || any(confidence <= 0) || any(confidence > 1)) stop("The confidence variable must be a numeric vector between 0 and 1")
	
	mcmcresults <- lapply(results, as.mcmc.list, add.mutate=TRUE)
	
	# Targets is a named list of monitor variables (which are automatically added to monitor) and the true value (if an array magically catch and rename to var[1] var[2] etc) which is compared to MCMC output - return distribution of medians and confidence% LCI/UCI and % within confidence limits - confidence may be a vector
	
	matchtargets <- function(sim) return(lapply(1:length(targets), function(x){
		# Allow a single variable name to match lots of sub-indexes (but not if it has square brackets already):
		if(!grepl("[",names(targets)[x],fixed=TRUE)){
			matching <- which(names(targets)[x] == gsub("\\[[[:digit:],]+\\]","",varnames(mcmcresults[[sim]])))
			ret <- unlist(rep(targets[[x]], length(matching)))
			names(ret) <- varnames(mcmcresults[[sim]])[matching]
		}else{
			ret <- targets[[x]]
			names(ret) <- names(targets)[x]
		}
		return(ret)
	}))
	
	matchedtargets <- matchtargets(1)
	
	if(!(all(sapply(1:length(results), function(sim) return(length(unlist(matchtargets(sim)))))==length(unlist(matchedtargets))))){
		stop("The specified named targets matched a different number of variables between simulations - ensure that any target arrays have equal dimensions between simulations", call.=FALSE)
	}	
	if(all(sapply(matchedtargets,length)==0)) stop(paste("None of the target variable(s) '", paste(names(targets),collapse="','"), "' were found in the JAGS output ('", paste(varnames(results[[1]]$mcmc),collapse="','"),"')",sep=""))
	if(any(sapply(matchedtargets,length)==0)) warning(paste("One or more of the target variable(s) '", paste(names(targets),collapse="','"), "' were not found in the JAGS output ('", paste(varnames(results[[1]]$mcmc),collapse="','"),"')",sep=""))
	
	targets <- unlist(matchedtargets)
	fullvars <- targets[alphabeticalvars(names(targets))]
	
	useindex <- which(varnames(mcmcresults[[1]]) %in% names(fullvars))
	usevars <- fullvars[names(fullvars) %in% varnames(mcmcresults[[1]])]
	
	usevars <- usevars[match(varnames(mcmcresults[[1]]),names(usevars))]
	usevars <- usevars[!is.na(usevars)]  # Need to exclude missing here - it's possible there are monitored variables that we don't care about

	if(any(is.na(usevars)) || !(all(names(usevars)==(varnames(mcmcresults[[1]])[useindex])))){
		stop(paste("Sorry - something has gone wrong with indexing the variables; I was expecting to see ", paste(names(usevars),collapse=","), " but saw ", paste((varnames(results[[1]]$mcmc)[useindex]),collapse=","), " - please submit a bug report!", sep=""))
	}
	
	simulations <- length(mcmcresults)

	sumtable <- vapply(mcmcresults,function(result){
		
		result <- combine.mcmc(result,collapse.chains=TRUE,vars=names(usevars))
		
		rtable <- matrix(nrow=length(useindex), ncol=4+(length(confidence)*4),dimnames=list(names(usevars), c("Target","Median","Mean",apply(expand.grid(c("Lower","Upper","Range","Within"),confidence*100,"%CI"),1,paste,collapse=""),"AutoCorr(Lag10)")))

		stochastic <- apply(result,2,var)!=0
		if(any(!stochastic)){
			if(all(!stochastic)) stop("All target variables are non-stochastic", call.=FALSE)
		}
		
		rtable[,1] <- usevars
		
		rtable[stochastic,2] <- apply(result[,stochastic,drop=FALSE],2,median)
		rtable[stochastic,3] <- apply(result[,stochastic,drop=FALSE],2,mean)
		for(c in 1:length(confidence)){
			h <- as.matrix(HPDinterval(result[,stochastic,drop=FALSE],prob=confidence[c]))
			r <- h[,2]-h[,1]
			w <- (usevars[stochastic] >= h[,1] & usevars[stochastic] <= h[,2])
			rtable[stochastic,((4*(c-1)) : ((4*c)-1))+4] <- c(h,r,w)
		}
		rtable[stochastic,ncol(rtable)] <- safe.autocorr.diag(result[,stochastic,drop=FALSE],lags=10)
		
		rtable[!stochastic,2:ncol(rtable)] <- NA
		
		return(rtable)
		
	}, matrix(0,nrow=length(useindex), ncol=4+(length(confidence)*4)))
	
	values <- apply(sumtable[,2,,drop=FALSE],1,function(x) return(sum(!is.na(x))))

#	Gives an inapporopriate warning with drop-1 type studies:	
	if(any(values==0) && runjags.getOption('summary.warning')) warning("One or more target variables is non-stochastic for all of the (successful) simulation results", call.=FALSE)
	
	meanable <- values>1
	if(!separatesingles)
		meanable[] <- TRUE
	
	if(any(meanable)){
		msum <- cbind(apply(sumtable[meanable,,,drop=FALSE],c(1,2),mean,na.rm=TRUE), Simulations=values[meanable])
		dimnames(msum) <- list(dimnames(msum)[[1]], c("Target","Av.Median","Av.Mean",apply(expand.grid(c("Av.Lower","Av.Upper","Av.Range","Prop.Within"),confidence*100,"%CI"),1,paste,collapse=""),"Av.AutoCorr(Lag10)","Simulations"))
	}else{
		msum <- NA	
	}
	indable <- values==1
	if(!separatesingles)
		indable[] <- FALSE
	if(any(indable)){
		isum <- t(apply(sumtable[indable,,,drop=FALSE],1,function(x) return(x[,which(!is.na(x[2,]))])))
	}else{
		isum <- NA	
	}
	
	# Re-add the non-stochastic targets here???
	
#	sumsum[,((4*(length(confidence)-1))+6)] <- sumsum[,((4*(length(confidence)-1))+6)]*100
	
	timetaken <- sapply(results, function(x) return(x$timetaken))
	sample <- sapply(results, function(x) return(x$sample * x$thin))
	burnin <- sapply(results, function(x) return(x$burnin * x$thin))
	
	return(list(means=msum, singles=isum, individual=sumtable, timetaken=timetaken, sample=sample, burnin=burnin))

}


run.JAGS.study <- run.jags.study