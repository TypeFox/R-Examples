#' @title Run or extend a user-specified Bayesian MCMC model in JAGS from within R
#' @name run.jags
#' @aliases run.jags run.JAGS extend.jags extend.JAGS
#' @export

#' @description
#' Runs or extends a user specified JAGS model from within R, returning an object of class \code{\link{runjags-class}}.  

#' @details
#' The run.jags function reads, compiles, and updates a JAGS model based on a model representation (plus data, monitors and initial values) input by the user.  The model can be contained in an external text file, or a character vector within R.  The autorun.jags function takes an existing \code{\link{runjags-class}} object and extends the simulation.  Running a JAGS model using these functions has two main advantages:
#' 
#' 1) The method used to call or extend the simulation can be changed simply using the method option.  The methods most likely to be used are 'interruptible' and 'rjags' which use one simulation per model, or 'parallel', 'bgparallel' and 'rjparallel' which run a separate simulation for each chain to speed up the model run.  For more details see below under the 'method' argument.  
#' 
#' 2) All information required to re-run the simulations is stored within the \code{\link{runjags-class}} object returned.  This complete representation can be written to a text file using \code{\link{write.jagsfile}}, then modified as necessary and re-run using only the file path as input.
#'
#' 3) Summary statistics for the returned simulations are automatically calculated and displayed using associated S3 methods intended to facilitate checking model convergence and run length.  Additional methods are available for plot functions, as well as conversion to and from MCMC and rjags objects.  See the help file for \code{\link{runjags-class}} for more details.  

#' @keywords models

#' @return
#' Usually an object of class 'runjags', or an object of class 'runjagsbginfo' for background methods (see \code{\link{runjags-class}}).

#' @seealso
#' \code{\link{results.jags}} to import completed simulations (or partially successful simulations) from saved JAGS files, \code{\link{runjags-class}} for details of available methods for the returned object, \code{\link{read.jagsfile}} for more details on the permitted format of the model file, \code{\link{write.jagsfile}} for a way to write an existing runjags object to file, and \code{\link{runjags.options}} for user options regarding warning messages etc.

#' @examples

#' \dontshow{
#' runjags.options(new.windows=FALSE)
#' }
#' # run a model to calculate the intercept and slope of the expression 
#' # y = m x + c, assuming normal observation errors for y:
#' 
#' # Simulate the data
#' X <- 1:100
#' Y <- rnorm(length(X), 2*X + 10, 1)
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
#' }"
#' 
#' # Data and initial values in a named list format, 
#' # with explicit control over the random number
#' # generator used for each chain (optional): 
#' data <- list(X=X, Y=Y, N=length(X))
#' inits1 <- list(m=1, c=1, precision=1,
#' .RNG.name="base::Super-Duper", .RNG.seed=1)
#' inits2 <- list(m=0.1, c=10, precision=1,
#' .RNG.name="base::Wichmann-Hill", .RNG.seed=2)
#' 
#' \dontrun{
#' # Run the model and produce plots 
#' results <- run.jags(model=model, monitor=c("m", "c", "precision"), 
#' data=data, n.chains=2, method="rjags", inits=list(inits1,inits2))
#' 
#' # Standard plots of the monitored variables:
#' plot(results)
#' 
#' # Look at the summary statistics:
#' print(results)
#' 
#' # Extract only the coefficient as an mcmc.list object:
#' library('coda')
#' coeff <- as.mcmc.list(results, vars="m")
#' }
#' 
#' 
#' # The same model but using embedded shortcuts to specify data, inits and monitors,
#' # and using parallel chains:
#' 
#' # Model in the JAGS format
#' 
#' model <- "model {
#' for(i in 1 : N){ #data# N
#' 	Y[i] ~ dnorm(true.y[i], precision) #data# Y
#' 	true.y[i] <- (m * X[i]) + c #data# X
#' }
#' m ~ dunif(-1000,1000) #inits# m
#' c ~ dunif(-1000,1000)
#' precision ~ dexp(1)
#' #monitor# m, c, precision
#' }"
#' 
#' # Simulate the data
#' X <- 1:100
#' Y <- rnorm(length(X), 2*X + 10, 1)
#' N <- length(X)
#' 
#' initfunction <- function(chain) return(switch(chain, 
#' 	"1"=list(m=-10), "2"=list(m=10)))
#' 
#' \dontrun{
#' # Run the 2 chains in parallel (allowing the run.jags function
#' # to control the number of parallel chains). We also use a
#' # mutate function to convert the precision to standard deviation:
#' results <- run.jags(model, n.chains=2, inits=initfunction,
#' method="parallel", mutate=list("prec2sd", vars="precision"))
#' 
#' # View the results using the standard print method:
#' results
#' 
#' # Look at some plots of the intercept and slope on a 3x3 grid:
#' plot(results, c('trace','histogram','ecdf','crosscorr','key'),
#' vars=c("m","^c"),layout=c(3,3))
#'
#' # Write the current model representation to file:
#' write.jagsfile(results, file='mymod.txt')
#' # And re-run the simulation from this point:
#' newresults <- run.jags('mymod.txt')
#' }
#' # Run the same model using 8 chains in parallel:
#' # distributed computing cluster:
#' \dontrun{
#'
#' # A list of 8 randomly generated starting values for m:
#' initlist <- replicate(8,list(m=runif(1,-20,20)),simplify=FALSE)
#' 
#' # Run the chains in parallel using JAGS (2 models 
#' # with 4 chains each):
#' results <- run.jags(model, n.chains=8, inits=initlist,
#' method="parallel", n.sims=2)
#'
#' # Set up a distributed computing cluster with 2 nodes:
#' library(parallel)
#' cl <- makeCluster(4)
#' 
#' # Run the chains in parallel rjags models (4 models 
#' # with 2 chains each) on this cluster:
#' results <- run.jags(model, n.chains=8, inits=initlist,
#' method="rjparallel", cl=cl)
#'
#' stopCluster(cl)
#'
#' # For more examples see the quick-start vignette:
#' vignette('quickjags', package='runjags')
#' 
#' # And for more details about possible methods see:
#' vignette('userguide', package='runjags')
#' }

#' @param model either a relative or absolute path to a textfile (including the file extension) containing a model in the JAGS language and possibly monitored variable names, data and/or initial values, or a character string of the same.  No default.  See \code{\link{read.jagsfile}} for more details.

#' @param monitor a character vector of the names of variables to monitor.  No default.  The special node names 'deviance', 'pd', 'popt', 'dic', 'ped' and 'full.pd' are used to monitor the deviance, mean pD, mean pOpt, DIC, PED and full distribution of sum(pD) respectively.  Note that these monitored nodes (with the exception of 'deviance') require multiple chains within the same simulation, and won't appear as variables in the summary statistics or plots (but see \code{\link{extract}} for a way of extracting these from the returned object).

#' @param data a named list, data frame, environment, character string in the R dump format (see \code{\link{dump.format}}), or a function (with no arguments) returning one of these types.  If the model text contains inline #data# comments, then this argument specifies the list, data frame or environment in which to search first for these variables (the global environment is always searched last).  If the model text does not contain #data# comments, then the full list or data frame (but not environment) is included as data.  If the data argument is a character string, then any #data# comments in the model are ignored (with a warning). The default specifies the parent environment of the function call.

#' @param n.chains the number of chains to use with the simulation.  More chains will improve the sensitivity of the convergence diagnostic, but will cause the simulation to run more slowly (although this may be improved by using a method such as 'parallel', 'rjparallel' or 'snow').  The minimum (and default) number of chains is 2.

#' @param inits either a character vector with length equal to the number of chains the model will be run using, or a list of named lists representing names and corresponding values of inits for each chain, or  a function with either 1 argument representing the chain or no arguments.   If a vector, each element of the vector must be a character string in the  R dump format representing the initial values for that chain, or NA.  If not all initialising variables are specified, the unspecified variables are taken deterministically from the mean or mode of the prior distribution by JAGS.  Values left as NA result in all initial values for that chain being taken from the prior distribution.  The special variables '.RNG.seed', '.RNG.name', and '.RNG.state' are allowed for explicit control over random number generators in JAGS.  If a function is provided, the data is available inside the function as a named list 'data' - this may be useful for setting initial values that depend on the data.  Default NA.  Note that the dimensions of any variables used for initial values must match the dimensions of the same parameter in the model - recycling is not performed.  If any elements of the initial values have deterministic values in the model, the corresponding elements must be defined as NA in the initial values.

#' @param runjags.object the model to be extended - the output of a run.jags (or autorun.jags or extend.jags etc) function, with class 'runjags'.  No default.

#' @param add.monitor a character vector of variables to add to the monitored variable list.  All previously monitored variables are automatically included - although see the 'drop.monitor' argument. Default no additional monitors.

#' @param drop.monitor a character vector of previously monitored variables to remove from the monitored variable list for the extended model. Default none.

#' @param drop.chain a numeric vector of chains to remove from the extended model. Default none.

#' @param combine a logical flag indicating if results from the new JAGS run should be combined with the previous chains.  Default TRUE if not adding or removing variables or chains, and FALSE otherwise.

#' @param burnin the number of burnin iterations, NOT including the adaptive iterations to use for the simulation.  Note that the default is 4000 plus 1000 adaptive iterations, with a total of 5000.

#' @param sample the total number of (additional) samples to take.  Default 10000 iterations.  If specified as 0, then the model will be created and returned without any MCMC samples (burnin and adapt will be ignored).  Note that a minimum of 100 samples is required to generate summary statistics.

#' @param adapt the number of adaptive iterations to use at the start of the simulation.  If the adaptive phase is not long enough, the sampling efficiency of the MCMC chains will be compromised.  If the model does not require adaptation (either because a compiled rjags model is already available or because the model contains no data) then this will be ignored, with a warning that the model is not in adaptive mode.  Default 1000 iterations.

#' @param datalist deprecated argument.

#' @param initlist deprecated argument.

#' @param noread.monitor an optional character vector of variables to monitor in JAGS and output to coda files, but that should not be read back into R.  This may be useful (in conjunction with keep.jags.files=TRUE) for looking at large numbers of variables a few at a time using the read.monitor argument to results.jags.  This argument is ignored for the rjags and rjparallel methods, and if keep.jags.files=FALSE.

#' @param jags the system call or path for activating JAGS.  Default uses the option given in \code{\link{runjags.options}}.

#' @param silent.jags option to suppress output of the JAGS simulations.  Default uses the option given in \code{\link{runjags.options}}.

#' @param modules a character vector of external modules to be loaded into JAGS, either as the module name on its own or as the module name and status separated by a space, for example 'glm on'.  

#' @param factories a character vector of factory modules to be loaded into JAGS.  Factories should be provided in the format '<facname> <factype> <status>' (where status is optional), for example: factories='mix::TemperedMix sampler on'.  You must also ensure that any required modules are also specified (in this case 'mix'). 

#' @param summarise should summary statistics be automatically calculated for the output chains?  Default TRUE (but see also ?runjags.options -> force.summary).

#' @param mutate either a function or a list with first element a function and remaining elements arguments to this function.  This can be used to add new variables to the posterior chains that are derived from the directly monitored variables in JAGS. This allows the variables to be summarised or extracted as part of the MCMC objects as if they had been calculated in JAGS, but without the computational or storage overheads associated with calculating them in JAGS directly.  The plot, summary and as.mcmc methods for runjags objects will automatically extract the mutated variables along with the directly monitored variables.  For an application to pairwise comparisons of different levels within fixed effects see \code{\link{contrasts.mcmc}}.

#' @param thin the thinning interval to be used in JAGS.  Increasing the thinning interval may reduce autocorrelation, and therefore reduce the number of samples required, but will increase the time required to run the simulation.  Using this option thinning is performed directly in JAGS, rather than on an existing MCMC object as with thin.sample. Default 1.

#' @param keep.jags.files option to keep the folder with files needed to call JAGS, rather than deleting it.  This allows the simulation results to be re-read using results.jags(path-to-folder), even from another R session, and may also be useful for attempting to bug fix models.   A character string can also provided, in which case this folder name  will be used instead of the default (existing folders will NOT be  over-written).  Default FALSE.  See also the \code{\link{cleanup.jags}} function.

#' @param tempdir option to use the temporary directory as specified by the system rather than creating files in the working directory.  If keep.jags.files=TRUE then the folder is copied to the working directory after the job has finished (with a unique folder name based on 'runjagsfiles').  Any files created in the temporary directory are removed when the function exits for any reason.  It is not possible to use a temporary directory with the background methods, so tempdir will be set to FALSE if not done so by the user (possibly with a warning  depending on the settings in \code{\link{runjags.options}}).  Default TRUE.

#' @param jags.refresh the refresh interval (in seconds) for monitoring JAGS output using the 'interactive' and 'parallel' methods (see the 'method' argument).  Longer refresh intervals will use slightly less processor time, but will make the simulation updates to be shown on the screen less frequently.  Reducing the refresh rate to every 10 or 30 seconds may be worthwhile for simulations taking several days to run.  Note that this will have no effect on the processor use of the simulations themselves.  Default 0.1 seconds.

#' @param batch.jags option to call JAGS in batch mode, rather than using input redirection.  On JAGS >= 3.0.0, this suppresses output of the status which may be useful in some situations.  Default TRUE if silent.jags is TRUE, or FALSE otherwise.

#' @param method the method with which to call JAGS; probably a character vector specifying one of 'rjags', 'simple', 'interruptible', 'parallel', 'rjparallel', 'background', 'bgparallel' or 'snow'. The 'rjags' and 'rjparallel' methods run JAGS using the rjags package, whereas other options do not require the rjags package and call JAGS as an external executable.  The advantage of the 'rjags' method is that the model will not need to be recompiled between successive calls to extend.jags, all other methods require a re-compilation (and adaptation if necessary) step at every call to extend.jags.  The 'background' and 'bgparallel' return a filename for the started simulation, which can be read using \code{\link{results.jags}}.  The 'parallel' and 'interruptible' methods for Windows require XP Professional, Vista or later (or any Unix-alike).  For more information refer to the userguide vignette.

#' @param method.options a deprecated argument currently permitted for backwards compatibility, but this will be removed from a future version of runjags.  Pass these arguments directly to run.jags or extend.jags.

#' @param ... summary parameters to be passed to \code{\link{add.summary}}, and/or additional options to control some methods including n.sims for parallel methods, cl for rjparallel and snow methods, remote.jags for snow, and by and progress.bar for the rjags method.
NULL


#' @rdname run.jags
run.jags <- function(model, monitor = NA, data=NA, n.chains=NA, inits = NA, burnin = 4000, sample = 10000, adapt=1000, noread.monitor=NULL, datalist=NA, initlist=NA, jags = runjags.getOption('jagspath'), silent.jags = runjags.getOption('silent.jags'), modules=runjags.getOption('modules'), factories=runjags.getOption('factories'), summarise = TRUE, mutate = NA, thin = 1, keep.jags.files = FALSE, tempdir=runjags.getOption('tempdir'), jags.refresh=0.1, batch.jags=silent.jags, method=runjags.getOption('method'), method.options=list(), ...){
	
	listwarn <- FALSE
	if(!identical(datalist, NA)){
		listwarn <- TRUE
		if(!identical(data, NA))
			stop('Arguments cannot be supplied for both datalist and data - please use the data argument only')
		data <- datalist
	}
	if(!identical(initlist, NA)){
		listwarn <- TRUE
		if(!identical(inits, NA))
			stop('Arguments cannot be supplied for both initlist and inits - please use the inits argument only')
		inits <- initlist
	}
	if(listwarn)
		warning('The datalist and initlist arguments are deprecated and will be removed from a future version of runjags - please use the data and inits arguments instead')
	
	# If data and inits are NA then grab the parent frame for first identification of variables:
	if(identical(data, NA))
		data <- parent.frame()
	if(is.null(inits) || identical(inits, NA))
		inits <- parent.frame()
	
	if(!identical(method.options, list()))
		swcat('Note:  the method.options argument is deprecated and will be removed from a future release of runjags\n')
	
	method <- getrunjagsmethod(method)
	obj <- setup.jagsfile(model=model, n.chains=n.chains, data=data, inits=inits, monitor=monitor, modules=modules, factories=factories, jags=jags, call.setup=TRUE, method=method, mutate=mutate)
	res <- extend.jags(runjags.object=obj, add.monitor=character(0), drop.monitor=character(0), drop.chain=numeric(0), combine=FALSE, burnin = burnin, sample = sample, adapt=adapt, noread.monitor=noread.monitor, jags = jags, silent.jags = silent.jags, summarise = summarise, thin = thin, keep.jags.files = keep.jags.files, tempdir=tempdir, jags.refresh=jags.refresh, batch.jags=batch.jags, method=method, method.options=method.options, ...)
	
	return(res)
	
}

#' @rdname run.jags
extend.jags <- function(runjags.object, add.monitor=character(0), drop.monitor=character(0), drop.chain=numeric(0), combine=length(c(add.monitor,drop.monitor,drop.chain))==0, burnin = 0, sample = 10000, adapt=1000, noread.monitor = NA, jags = NA, silent.jags = NA, summarise = sample >= 100, thin = NA, keep.jags.files = FALSE, tempdir=runjags.getOption('tempdir'), jags.refresh=NA, batch.jags=silent.jags, method=NA, method.options=NA, ...){
	
	runjags.object <- checkvalidrunjagsobject(runjags.object)
	
	# We get unhelpful error messages otherwise:
	if(identical(NA, method))
		method <- runjags.object$method	
	if(identical(NA, method.options) || identical(list(), method.options)){
		method.options <- runjags.object$method.options
	}else{
		# Needed to re-create default options if we are given a partial list - allow nsims for backwards compatibility
		nmo <- method.options
		method.options <- getdefaultmethodoptions()
		nmonam <- names(nmo)
		if(any(nmonam=="nsims") && !any(nmonam=="n.sims")){
			nmonam[nmonam=="nsims"] <- "n.sims"
			names(nmo) <- nmonam
		}
		for(n in names(nmo))
			method.options[n] <- nmo[n]
	}

	if(identical(NA, summarise))
		stop('NA value not allowed for summarise argument')
	
	if(identical(NA, noread.monitor))
		noread.monitor <- runjags.object$noread.monitor
	if(identical(NA, jags))
		jags <- runjags.object$method.options$jags
	if(identical(NA, silent.jags))
		silent.jags <- runjags.object$method.options$silent.jags
	if(identical(NA, thin))
		thin <- runjags.object$thin
	if(identical(NA, jags.refresh))
		jags.refresh <- runjags.object$method.options$jags.refresh
	if(identical(NA, batch.jags))
		batch.jags <- runjags.object$method.options$batch.jags
	
	if(jags.refresh > 60)
		warning('You have set jags.refresh to a very low frequency - this will not appreciably reduce the computational effort required and may unnecessarily delay loading the results')
		
	# We may be passed some unevaluated function arguments from parent functions using getargs so evaluate everything here:
	argnames <- names(formals(extend.jags))
	argnames <- argnames[argnames!='...']
	for(i in 1:length(argnames)){
		success <- try(assign(argnames[i], eval(get(argnames[i]))), silent=TRUE)		
		if(inherits(success, 'try-error')){
			stop(paste("object '", strsplit(as.character(success),split="'",fixed=TRUE)[[1]][2], "' not found", sep=""), call.=FALSE)
		}
	}
	
	# Check summary arguments before running model:
	summaryargs <- getsummaryargs(runjags.object$summary.pars, 'extend.jags', silent.jags=silent.jags, ignore=c("n.sims","cl","remote.jags","rjags","by","progress.bar"), ...)
	passed <- list(...)
	
	# Have to stop, otherwise it will crash later - from autorun that now doesn't check ...
	if(any(names(passed)=='adaptive')){
		warning('The adaptive argument is deprecated and has been removed - please use the adapt argument instead')
		passed$adaptive <- NULL
	}
	
	# Combine code won't work if there is nothing to combine, so pretend its new:
	if(runjags.object$sample==0)
		combine <- FALSE
	if(is.na(runjags.object$method)){
#		warning("Cannot combine old and new MCMC simulations when given an empty initialised runjags object - setting combine to FALSE")		
#		combine <- FALSE
		method <- if('rjags' %in% .packages(TRUE)) 'rjags' else 'interruptible'
	}
	method <- getrunjagsmethod(method)	
	n.chains <- length(runjags.object$end.state)	
	
	########################## Set method options (and what will become method options) here:
	##### Defulat method options are set by setup.jags	

	# Currently ignoring if these options are available but eventually warn if the method doesn't match:
	#warning("The supplied value for method.options is ignored when using an inbuilt method (except for 'n.sims' which can be provided for parallel methods, 'cl' and 'remote.jags' which can be provided for the snow method, and 'by' and 'progress.bar' which can be supplied for the rjags method)")
	# parallel:  c('n.sims')
	# snow:  c('cl', 'remote.jags', 'n.sims')
	# all not rjags/rjparallel:  batch.jags, (also kind of jags, jags.refresh, tempdir and keep.jags.files - change in future?)
	# rjags:  c('by', 'progress.bar')
	# rjparallel:  c('by', 'progress.bar', 'cl', 'n.sims')

	if('cl' %in% names(passed))
		method.options$cl <- passed$cl
	if(!identical(method.options$cl, NA) && !any(grepl('cluster',class(method.options$cl))))
		stop('Unrecognised argument for "cl" - this should be a cluster')

	if('n.sims' %in% names(passed))
		method.options$n.sims <- passed$n.sims
	if(!identical(method.options$n.sims, NA) && (!is.numeric(method.options$n.sims) || length(method.options$n.sims)!=1))
		stop('Unrecognised argument for "n.sims" - this should be a numeric')
		
	if('remote.jags' %in% names(passed))
		method.options$remote.jags <- passed$remote.jags

	if('by' %in% names(passed))
		method.options$by <- passed$by

	if('progress.bar' %in% names(passed))
		method.options$progress.bar <- passed$progress.bar
	
	passed <- passed[! names(passed) %in% c("n.sims","cl","remote.jags","rjags","by","progress.bar")]
	
	# Some of these are used so will have to be visible, but they really should be part of this:
	method.options$jags <- jags
	method.options$silent.jags <- silent.jags
	method.options$jags.refresh <- jags.refresh
	method.options$batch.jags <- batch.jags
	########################## 
	
	# Check modules and factories
	runjags.object$factories <- checkmodfact(runjags.object$factories, 'factory')
	runjags.object$modules <- checkmodfact(runjags.object$modules, 'module')
	
	if(length(keep.jags.files)!=1) stop("keep.jags.files must be either TRUE, FALSE or a character string specifying a folder to save results to")
	if(class(tempdir)!="logical") stop("tempdir must be either TRUE or FALSE")
	if(class(keep.jags.files)=="character"){
		tempdir <- FALSE
		directory <- keep.jags.files
		keep.jags.files <- TRUE
	}else{
		directory <- NA
	}
	
	starttime <- Sys.time()
		
	# Adjust monitors and chains:

	if(combine && runjags.object$thin != thin)
		stop("Cannot combine old and new MCMC simulations when changing the thin interval")
	if(length(c(add.monitor,drop.monitor,drop.chain))!=0 && combine) stop("Cannot combine old and new MCMC simulations when adding or removing monitors or chains")
	
	monitor <- runjags.object$monitor
	modules <- runjags.object$modules
		
	if(sample==0){
#		if((adapt>0 && formals(extend.jags)$adapt!=adapt) || (burnin>0 && formals(extend.jags)$burnin!=burnin))
#			warning('Cannot adapt or burn in the model with sample=0 - setting adapt=0 and burnin=0')
		adapt <- 0
		burnin <- 0
		summarise <- FALSE
		# DO NOT set the method to simple here - it breaks autorun.jags ... if the end user wants this simple, they specify method='simple'
		# method <- 'simple'
	}
		
	if(sample < 0) stop("A positive integer must be supplied for 'sample'")
	if(summarise && ((sample + combine*runjags.object$sample) < 100) && !runjags.getOption('force.summary')){
		if(runjags.getOption('summary.warning'))
			warning("Cannot produce meaningful summary statistics with less than 100 samples", call.=FALSE)
		summarise <- FALSE
	} 

#	if(!summarise && length(passed)!=0 && runjags.getOption('summary.warning'))
#		warning('Options to add.summary are ignored when summarise=FALSE', call.=FALSE)
	
	if(n.chains > 1 && all(runjags.object$end.state=="")){
		if(runjags.getOption('inits.warning')) warning("No initial values were provided - JAGS will use the same initial values for all chains", call.=FALSE)
	}
	
	# In case this has been changed by the user (ignore potential for drop monitor to indicate a noread.monitor):
	runjags.object$noread.monitor <- checkvalidmonitorname(noread.monitor)
	
	if(length(drop.monitor)>0){
		drop.monitor <- checkvalidmonitorname(drop.monitor)
		remove <- which(monitor %in% drop.monitor)
		if(length(remove)!=length(drop.monitor)) warning(paste("Specified variables ", paste(drop.monitor[!(drop.monitor %in% monitor)], collapse=", "), " were not being monitored so could not be dropped", sep=""))
		if(length(remove)>0) monitor <- monitor[-remove]
	}	
	
	if(length(add.monitor)>0){
		add.monitor <- checkvalidmonitorname(add.monitor)
		monitor <- c(monitor, add.monitor)

		monitor[monitor=="DIC"] <- "dic"
		monitor[monitor=="PED"] <- "ped"
		monitor[monitor=="pD"] <- "pd"
		monitor[monitor=="full.pD"] <- "full.pd"
		monitor[monitor=="pOpt"] <- "popt"	

		if(any(c("pd","full.pd")%in%monitor)){
			monitor <- c(monitor, "dic")
		}
		if(any(monitor=="popt")){
			monitor <- c(monitor, "ped")
		}
	
		monitor[monitor=="pD.i"] <- "pd.i"
		monitor <- na.omit(monitor[monitor!=""])
		monitor <- unique(monitor)
	}
	if(any(monitor=='pd.i')){
		warning('Ignoring deprecated pd.i monitor', call.=FALSE)
		monitor <- monitor[monitor!='pd.i']
	}
			
	
	if(any(noread.monitor %in% monitor)){
		warning('Duplicate variables detected in monitor and noread.monitor - names specified to monitor take precedence')
		noread.monitor <- noread.monitor[!noread.monitor%in%monitor]
	}
	noread.monitor <- checkvalidmonitorname(noread.monitor)		
		
	if(any(c("deviance","ped","dic") %in% monitor) && as.character(runjags.object$data)==""){
		warning("Unable to monitor deviance with no data - removing any deviance/DIC etc monitors")
		monitor <- monitor[!monitor%in%c("deviance","pd","full.pd","popt","dic","ped")]
	}
	
 	modules <- checkmodfact(modules, 'module')
	forcerecompile <- FALSE
	if(any(c("popt", "pd", "full.pd", "deviance", "dic", "ped") %in% monitor)){
		if(identical(modules,'') || !any('dic'==sapply(modules,function(x) return(x[1]))))
			forcerecompile <- TRUE
		modules <- c(modules, list(c("dic","TRUE")))
	}
	if(any(grepl("lecuyer::RngStream",runjags.object$end.state))){
		if(identical(modules,'') || !any(unlist(modules)=='lecuyer'))
			forcerecompile <- TRUE
		modules <- c(modules, list(c("lecuyer","TRUE")))
	}
 	modules <- checkmodfact(modules, 'module')
	
	monitor[monitor==""] <- NA
	if(class(monitor)!="character" | all(is.na(monitor))){
		stop("Monitored variable(s) must be provided in the form of a character vector")
	}
	monitor <- unique(na.omit(monitor))

	if(length(drop.chain)>0){

		drop.chain <- unique(round(drop.chain))
		if(any(drop.chain>length(runjags.object$end.state)) | any(drop.chain < 1)) stop("Specified value(s) to drop.chain are invalid - please specify which chain(s) to drop by the chain number(s)")
		if(length(drop.chain)==length(runjags.object$end.state)) stop("Specified value(s) to drop.chains are invalid - it is not possible to drop all chains")

		# Will have to re-compile rjags object:
		method.options <- method.options[names(method.options)!='rjags']
		runjags.object$end.state <- runjags.object$end.state[-drop.chain]		
	}

	inits <- runjags.object$end.state	
	n.chains <- length(inits)
	if(n.chains<1) stop("Number of chains must be greater than 0")
		
		
	if(as.integer(thin)!=thin | thin < 1) stop("The value supplied for thin must be a positive integer")
	
	sample <- sample
	burnin <- burnin
	
	
	###########
	### Setup rjags stuff:
	# Check to see if this is using an rjags method, and if it is get the method.options$rjags stuff set up:	
	if(method%in%runjagsprivate$parallelmethod && any(c('ped','dic','pd','full.pd','popt')%in%monitor))
		stop("The DIC, PED, pD, full.pD and pOpt cannot be assessed when using parallel or separate chains")
	
	if(length(noread.monitor)>0 && method%in%runjagsprivate$rjagsmethod){
		warning("The noread.monitor argument cannot be used with the 'rjags' or 'rjparallel' methods and was ignored")
		noread.monitor <- NULL
	}
	if(keep.jags.files && method%in%runjagsprivate$rjagsmethod){
		warning("Unable to keep JAGS files when using the 'rjags' or 'rjparallel' methods - switching to the 'interruptible' or 'parallel' method")
		method <- if(method=="rjags") "interruptible" else "parallel"
	}
	if(tempdir && method%in%runjagsprivate$bgmethod){
		if(runjags.getOption('tempdir.warning')) warning("Working in a temporary directory is not advisable for background methods - creating JAGS files in the current working directory")
		tempdir <- FALSE
	}
		
	# Removed because BANOVA specifies adapt=0 - it does make sense under some circumstances, and a warning is given if not adapted
	# if((adapt+burnin+sample)>0 && adapt < 1) warning('An adaptation phase of >0 iterations must be specified unless extending an already compiled model using the rjags method', call.=FALSE)	
	
	if(method%in%runjagsprivate$rjagsmethod){
		
    	# RNG in >4 rjparallel chains should be lecuyer - if the lecuyer module is loaded when the rjags object compiles (unless it already has RNGnames) it will be
		if(n.chains > 4 && method=='rjparallel'){
						
			norng <- sum(grepl('.RNG.name', inits, fixed=TRUE))
			if(norng!=n.chains && norng!=0) stop('Attempting to use parallel chains with some (but not all) .RNG.name values specified - make sure you specify a .RNG.name for each chain (or no chains) and try again')
			# If we have RNG names make sure they are all lecuyer:
			if(norng!=0){
				if(sum(grepl('lecuyer::RngStream', inits, fixed=TRUE)) != norng){
					if(runjags.getOption('rng.warning')) warning("You attempted to start >4 parallel chains with a PRNG other than that in the lecuyer module - this is not recommended.  The RNG.names have been modified to use the lecuyer module.")
					inits <- lapply(inits,function(x){
						x <- list.format(x)
						x <- x[!names(x)%in%c(".RNG.state",".RNG.name")]
						return(dump.format(x))
					})
					class(inits) <- "runjagsinits"
					runjags.object$end.state <- inits
					method.options <- method.options[names(method.options)!='rjags']
					# The inits now don't have .RNG.name and we have loaded the lecuyer module - so when the model re-compiles (as it will be forced to), the new inits will be set as lecuyer
				}
			}else{
        		# Otherwise make sure the lecuyer module is loaded:
				modules <- c(modules, list(c("lecuyer","TRUE")))
			}
		}
		
		runjags.object$modules <- checkmodfact(modules, 'module')
		runjags.object$method.options <- method.options

		if('rjags' %in% names(method.options)){
      		# Checks compiled:
			checkcompiled <- try(stats::coef(runjags.object$method.options$rjags),silent=TRUE)
			
			# If repeatable methods force a recompile:
			if(runjags.getOption('repeatable.methods') || forcerecompile || class(checkcompiled)=='try-error'){
				runjags.object$method.options <- runjags.object$method.options[names(runjags.object$method.options)!='rjags']
				method.options$rjags <- as.jags(runjags.object, adapt=0, quiet=silent.jags)
			}
			#else{
				# We could check that modules and factories we need are already loaded here and load them if not - I don't actually think recompiling is necessary (but test)
	#			forcerecompile <- FALSE
	#			loadedmods <- list.modules()
	#			for(m in modules){
	#				if(m!=''){
	#					if()
	#				}
	#			}
				
	#    It is a bad idea to override adapt here - if it isn't required, it just won't be done:
				# ONLY if using rjags AND NOT repeatable methods - otherwise the model will be recompiled:
#				if(method=='rjags')
#         			adapt <- 0

			#}

		}else{
			# Sets up and compiles:
			method.options <- c(method.options, list(rjags=as.jags(runjags.object, adapt=0, quiet=silent.jags)))
		}
				
	}else{
		method.options <- method.options[names(method.options)!="rjags"]
	
		# None of this will ever be needed for rjags methods - rjparallel separate chains are handled above as it needs to be done before model compilation:
		if(method%in%runjagsprivate$parallelmethod){

			norng <- sum(grepl('.RNG.name', inits, fixed=TRUE))

			if(norng!=n.chains){
				if(norng!=0) 
					stop('Attempting to use parallel chains with some (but not all) .RNG.name values specified - make sure you specify a .RNG.name for each chain and try again', call.=FALSE)
				if(runjags.getOption('rng.warning')) 
					warning('You attempted to start parallel chains without setting different PRNG for each chain, which is not recommended.  Different .RNG.name values have been added to each set of initial values.', call.=FALSE)
		
				# Different behaviours for nchains <= 4 and > 5:

				if(n.chains <= 4){
					rngname <- paste('\".RNG.name\" <- \"', rep(c('base::Wichmann-Hill', 'base::Marsaglia-Multicarry', 'base::Super-Duper', 'base::Mersenne-Twister'), ceiling(n.chains/4)), '\"\n',sep='')
				}else{
					rjagsloaded <- 'rjags' %in% .packages()
					if(!loadandcheckrjags(FALSE))
						stop("The rjags package is required to generate initial values for more than 4 parallel chains")
					
					success <- try(rjags::load.module("lecuyer"))
					if(inherits(success, 'try-error')) stop("Failed to load the lecuyer module - ensure that the latest version of JAGS and rjags is installed")
				
					rngname <- sapply(rjags::parallel.seeds("lecuyer::RngStream", n.chains), dump.format)		
					if(!any(sapply(modules,function(x) return(x[1]))=='lecuyer'))
						modules <- c(modules, list(c("lecuyer","TRUE")))
				
					modules <- checkmodfact(modules, 'module')
				
					if(!rjagsloaded) unloadNamespace('rjags')		
				}
					
				for(i in 1:n.chains){
					if(is.na(inits[i])) inits[i] <- rngname[i] else inits[i] <- paste(inits[i], '\n', rngname[i], sep='')
				}		

			}
		
		}
	}
		
	###########
	###########
	
	if(sample==0 && (!combine || runjags.object$sample==0)){
		# Return a runjags object with 0 iterations:
		timetaken <- (difftime(Sys.time(), starttime) + runjags.object$timetaken)
		
		blankmcmc <- as.mcmc.list(lapply(1:n.chains, function(x) return(as.mcmc(matrix(NA,ncol=1,nrow=0)))))			
		
		# Destroy the compiled rjags model (if any) - it hasn't been adapted:
		if((runjags.object$sample+runjags.object$burnin)==0)
			method.options <- method.options[names(method.options)!="rjags"]
		
		combinedoutput <- makerunjagsobject(list(mcmc=blankmcmc, deviance.table=NA, deviance.sum=NA, pd=NA, end.state=inits, samplers=NA, end.time=Sys.time()), summarise=FALSE, summaryargs=summaryargs, burnin=runjags.object$burnin+runjags.object$sample, sample=0, thin=thin, model=runjags.object$model, data=runjags.object$data, monitor=monitor, noread.monitor=noread.monitor, modules=modules, factories=runjags.object$factories, response=runjags.object$response, residual=runjags.object$residual, fitted=runjags.object$fitted, method=method, method.options=method.options, timetaken=(difftime(Sys.time(), starttime) + timetaken), silent=TRUE)
		
		return(combinedoutput)
	}
	
	# Wrapper to catch sample=0, in which case we are probably just using this as a hack to calculate sumamry statistics:
	if(sample > 0){
		# Call function to run simulation and return MCMC list and pd object:
		startinfo <- runjags.start(model=runjags.object$model, monitor=c(monitor, noread.monitor), data=runjags.object$data, inits=inits, modules=modules, factories=runjags.object$factories, burnin = burnin, sample = sample*thin, adapt=adapt, thin = thin, tempdir=tempdir, dirname=directory, method=method, method.options=method.options)
		
		# If using rjparallel method destroy the compiled rjags model as we update different compiled objects:
		if(method=='rjparallel' || runjags.getOption('repeatable.methods'))
			method.options <- method.options[names(method.options)!="rjags"]
		
		# Used when there is an error, or when the simulations are background:
		dumpmodelinfo <- function(save.directory, returnbg){
		
			# Make a background.runjags.object object and save it - to be used by user in case of emergencies
			background.runjags.object <- runjags.object
			if(!combine) background.runjags.object$mcmc <- NA
			background.runjags.object$pd <- NA
			background.runjags.object$deviance.table <- NA
			background.runjags.object$deviance.sum <- NA
			background.runjags.object$end.state <- NA
			background.runjags.object$thin <- thin
			background.runjags.object$HPD <- NA
			background.runjags.object$hpd <- NA
			background.runjags.object$mcse <- NA
			background.runjags.object$psrf <- NA
			background.runjags.object$autocorr <- NA
			background.runjags.object$crosscorr <- NA
			background.runjags.object$truestochastic <- NA
			background.runjags.object$semistochastic <- NA
			background.runjags.object$nonstochastic <- NA
			background.runjags.object$dic <- NA
			background.runjags.object$trace <- NA
			background.runjags.object$density <- NA
			background.runjags.object$thin <- NA
		
			background.runjags.object$monitor <- monitor
			background.runjags.object$noread.monitor <- noread.monitor
			background.runjags.object$method <- method
			background.runjags.object$method.options <- method.options
			
			background.runjags.object$combine <- combine
			background.runjags.object$summarise <- summarise
			background.runjags.object$silent.jags <- silent.jags
			background.runjags.object$oldburnin <- runjags.object$burnin
			background.runjags.object$oldsample <- runjags.object$sample
			background.runjags.object$adapt <- adapt
			background.runjags.object$burnin <- burnin
			background.runjags.object$sample <- sample
			background.runjags.object$thin <- thin
			background.runjags.object$oldthin <- runjags.object$thin
			background.runjags.object$keep.jags.files <- keep.jags.files
			background.runjags.object$tempdir <- tempdir
			background.runjags.object$inits <- inits
			background.runjags.object$modules <- modules
			
			background.runjags.object$summaryargs <- summaryargs
		
			background.runjags.object$startedon <- starttime
			background.runjags.object$runjags.version <- c(runjagsprivate$runjagsversion, R.Version()$version.string, .Platform$OS.type, .Platform$GUI, .Platform$pkgType, format(Sys.time()))
		
			background.runjags.object <- c(background.runjags.object, startinfo)
		
			class(background.runjags.object) <- "runjagsbginfo"
		
			save(background.runjags.object, file=file.path(save.directory,"jagsinfo.Rsave"))
			
			# Return an object of class suspendedrunjags with all these options inside it, and the part of the old runjags object that we still need (eg for joining etc) also inside along with startinfo itself
			if(returnbg) return(background.runjags.object) else return(NULL)
			
		}
		
		if(!startinfo$complete){
			background.runjags.object <- dumpmodelinfo(startinfo$directory, TRUE)
			runjagsprivate$simfolders <- c(runjagsprivate$simfolders, startinfo$directory)												
			return(background.runjags.object)
		}
		
		
		# If using rjags or rjparallel just get the MCMC stuff from there:	
		if(all(c("mcmc","end.state") %in% names(startinfo))){
			if(!"pd" %in% names(startinfo)){
				startinfo$pd <- NA
			}
			if(any(monitor=="full.pd") && identical(startinfo$pd, NA)){
				warning("The full.pD was not returned as expected; removing the full.pD monitor")
				monitor <- monitor[! monitor %in% c("full.pd")]
			}
			newoutput <- list(mcmc=startinfo$mcmc, deviance.table=startinfo$deviance.table, deviance.sum=startinfo$deviance.sum, pd=startinfo$pd, end.state=startinfo$end.state, samplers=startinfo$samplers)
			
		}else{
			# Otherwise we have to call runjags.readin and then deal with copying files etc:
			
			if(length(noread.monitor)>0){
				read.monitor <- monitor[!monitor%in%c('pd','full.pd','popt','dic')]
			}else{
				read.monitor <- NA
			} 
			
			
			allok <- FALSE
			# This will be called if runjags.readin fails (allok=FALSE) or at the end of the function (allok=TRUE):
			on.exit({
				
				new.directory <- startinfo$directory
				if(keep.jags.files && tempdir){
					new.directory <- if(class(method)=="list" && method$method=='xgrid') new_unique(method$jobname, touch=TRUE, type='folder') else new_unique('runjagsfiles', touch=TRUE, type='folder')			
					if(new.directory=="Directory not writable"){
						warning("JAGS files could not be copied to the working directory as it is not writable")
					}else{
						# file.rename may not work on all platforms for directories so use file.copy instead:
						file.copy(from=file.path(startinfo$directory, list.files(startinfo$directory)), to=new.directory, recursive=TRUE)
						swcat("JAGS files were saved to the '", new.directory, "' folder in your current working directory\n", sep="")
					}
				}
        
				if(keep.jags.files){
					if(!allok)
					  swcat('Note: Either one or more simulation(s) failed, or there was an error in processing the results.  You may be able to retrieve any successful simulations using:\nresults.jags("', new.directory, '", recover.chains=TRUE)\nSee the help file for that function for possible options.\n', sep='')
          
					runjagsprivate$simfolders <- c(runjagsprivate$simfolders, new.directory)									
				}else{
					if(!allok && runjags.getOption('keep.crashed.files') && new.directory!="Directory not writable"){
					  swcat('Note: Either one or more simulation(s) failed, or there was an error in processing the results.  You may be able to retrieve any successful simulations using:\nresults.jags("', new.directory, '", recover.chains=TRUE)\nSee the help file for that function for possible options.\n', sep='')					  
						swcat('To remove failed simulation folders use cleanup.jags() - this will be run automatically when the runjags package is unloaded\n')
						runjagsprivate$failedsimfolders <- c(runjagsprivate$failedsimfolders, new.directory)
						keep.jags.files <- TRUE
					}
				}

				# If we either want to delete files (and the simulation didn't crash), or we have copied from the temp dir, delete the sim folder:
				if(!keep.jags.files || (allok && tempdir)) unlink(startinfo$directory, recursive = TRUE)
					
				# Now create the JAGS model info in the save folder if keeping files (the return value is ignored):
				if(keep.jags.files) dumpmodelinfo(new.directory, returnbg=FALSE)
				
			})
			
			newoutput <- runjags.readin(directory=startinfo$directory, silent.jags=silent.jags, target.adapt=adapt, target.burnin=burnin, target.iters=sample, n.chains=length(inits), monitor=monitor, method=method, method.options=method.options, suspended=FALSE, read.monitor=read.monitor, sub.samples=FALSE, sub.chains=if(runjags.getOption('partial.import') && !combine) TRUE else FALSE, force.read=length(noread.monitor)>0)

			# The time between starting and the sims finishing:
#			timetaken <- (difftime(newoutput$end.time, starttime) + runjags.object$timetaken)
			
			# unfinished will be NULL or TRUE - but should always be NULL for extend.jags
			if(identical(newoutput$unfinished,TRUE)){
				stop('An unexpected error occured - the JAGS simulations appear to not be finished.  Please file a bug report to the package author.')
			}

			allok <- TRUE			

			if('errormessage'%in%names(newoutput))
				stop(newoutput$errormessage,call.=FALSE)
			
		}
		
		end.state <- newoutput$end.state
		burnin <- runjags.object$burnin+(runjags.object$thin*runjags.object$sample)+burnin+adapt
		
		iternames <- seq((burnin+1), (burnin+(sample*thin))-(thin-1), length.out=sample)
		currentdn <- dimnames(newoutput$mcmc[[1]]) 
    
		# rjags seems to count adapt before burnin, but separate JAGS doesn't.  Also, rjparallel starts from iteration 1.
		#if(!all(currentdn[[1]]==iternames) && !all(currentdn[[1]]==(iternames-adapt))){
		#  warning('The iteration labels for the returned MCMC object were not as expected and have been renamed')
		#}
		
		for(i in 1:length(newoutput$mcmc)){
			if(dim(newoutput$mcmc[[i]])[1]!=sample) stop('An unexpected error occured when renaming and combining the MCMC chains')
			newoutput$mcmc[[i]] <- mcmc(newoutput$mcmc[[i]], start=burnin+1, thin=thin)
			dimnames(newoutput$mcmc[[i]]) <- list(iternames, currentdn[[2]])
		}
		if(class(newoutput$pd)=="mcmc" && !is.na(newoutput$pd)){
			dimnames(newoutput$pd) <- list(iternames, dimnames(newoutput$pd)[[2]])
			newoutput$pd <- mcmc(newoutput$pd, start=burnin+1, thin=thin)		
		}	
		
		# Combine runjags objects if necessary:
		if(combine){
			s <- try(combpd <- if('full.pd'%in%monitor) combine.mcmc(list(runjags.object$pd, newoutput$pd), collapse.chains=FALSE) else NA)
			if(class(s)=='try-error'){
				warning('An unexpected error occured while trying to combine the full.pD from the old and new simulations - it has been removed', call.=FALSE)
				monitor <- monitor[monitor!='full.pd']
				combpd <- NA
			}
			combinedoutput <- list(mcmc=combine.mcmc(list(runjags.object$mcmc, newoutput$mcmc), collapse.chains=FALSE), deviance.table=weightedaverage(runjags.object$deviance.table, newoutput$deviance.table, niter(runjags.object$mcmc), niter(newoutput$mcmc)), deviance.sum=weightedaverage(runjags.object$deviance.sum, newoutput$deviance.sum, niter(runjags.object$mcmc), niter(newoutput$mcmc)), pd=combpd, end.state=newoutput$end.state, samplers=newoutput$samplers)

		}else{
			combinedoutput <- list(mcmc=newoutput$mcmc, deviance.table=newoutput$deviance.table, deviance.sum=newoutput$deviance.sum, pd=if('full.pd'%in%monitor) newoutput$pd else NA, end.state=newoutput$end.state, samplers=newoutput$samplers)
		}
		
		# Save some RAM:
		rm(newoutput)

	}else{
		combinedoutput <- list(mcmc=runjags.object$mcmc, deviance.table=runjags.object$deviance.table, deviance.sum=runjags.object$deviance.sum, pd=runjags.object$pd, end.state=runjags.object$end.state, samplers=runjags.object$samplers)			
	}
	
	if(combine) burnin <- runjags.object$burnin
	
	# For submit and stop use the file modification date returned by runjags.readin
	timetaken <- (difftime(Sys.time(), starttime) + runjags.object$timetaken)
	
	combinedoutput <- makerunjagsobject(combinedoutput, summarise=summarise, summaryargs=summaryargs, burnin=burnin, sample=niter(combinedoutput$mcmc), thin=thin, model=runjags.object$model, data=runjags.object$data, monitor=monitor, noread.monitor=noread.monitor, modules=modules, factories=runjags.object$factories, response=runjags.object$response, residual=runjags.object$residual, fitted=runjags.object$fitted, method=method, method.options=method.options, timetaken=timetaken)
	
	if(sample > 0)
		swcat("Finished running the simulation\n")
	
	stopifnot(class(combinedoutput$end.state)=='runjagsinits')
	
	
#	if(identical(combinedoutput$samplers, NA) && runjags.getOption('nodata.warning'))
#		warning('JAGS is reporting that no samplers were used - ensure that any data has been passed to JAGS correctly')
#  This can happen when just forward sampling from the prior
	
	return(combinedoutput)
}


run.JAGS <- run.jags
extend.JAGS <- extend.jags
