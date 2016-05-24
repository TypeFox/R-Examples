#' @title Run or extend a user-specified Bayesian MCMC model in JAGS with automatically calculated run-length and convergence diagnostics
#' @name autorun.jags
#' @aliases autorun.jags autorun.JAGS autoextend.jags autoextend.JAGS
#' @export

#' @description
#' Runs or extends a user specified JAGS model from  within R, returning an object of class \code{\link{runjags-class}}.  The model is automatically assessed for convergence and adequate sample size before being returned.

#' @details
#' The autorun.jags function reads, compiles, and updates a JAGS model based on a model representation (plus data, monitors and initial values) input by the user.  The autoextend.jags function takes an existing \code{\link{runjags-class}} object and extends the simulation as required.  Chain convergence over the first run of the simulation is assessed using Gelman and Rubin's convergence diagnostic.  If necessary, the simulation is extended to improve chain convergence (up to a user-specified maximum time limit), before the required sample size of the Markov chain is calculated using Raftery and Lewis's diagnostic.  The simulation is extended to the required sample size dependant on autocorrelation and the number of chains. Note that automated convergence diagnostics are not perfect, and should not be considered as a replacement for manually assessing convergence and Monte Carlo error using the results returned.  For more complex models, the use of \code{\link{run.jags}} directly with manual assessment of necessary run length may be preferable.
#'
#' For autoextend.jags, any arguments with a default of NA are taken from the runjags object passed to the function.

#' @keywords models

#' @return
#' An object of class 'runjags' (see \code{\link{runjags-class}} for available methods).

#' @seealso
#' \code{\link{run.jags}} for fixed run length models, \code{\link{read.winbugs}} for details of model specification options, \code{\link{read.jagsfile}} and \code{\link{summary.runjags}} for details on available methods for the returned models, and \code{\link{run.jags.study}} for examples of simulation studies using automated model control provided by autorun.jags

#' @examples
#' # Run a model to calculate the intercept and slope of the expression 
#' # y = m x + c, assuming normal observation errors for y:
#' 
#' # Simulate the data
#' N <- 100
#' X <- 1:N
#' Y <- rnorm(N, 2*X + 10, 1)
#' 
#' # Model in the JAGS format
#' model <- "model {
#' for(i in 1 : N){
#' 	Y[i] ~ dnorm(true.y[i], precision)
#' 	true.y[i] <- m * X[i] + c
#' }
#' m ~ dunif(-1000,1000)
#' c ~ dunif(-1000,1000)
#' precision ~ dexp(1)
#' 
#' #data# N, X, Y
#' #inits# m, c, precision
#' }"
#' 
#' # Initial values to be used:
#' m <- list(-10, 10)
#' c <- list(-10, 10)
#' precision <- list(0.1, 10)
#' \dontrun{
#' # Run the model using rjags with a 5 minute timeout:
#' results <- autorun.jags(model=model, max.time="5m", 
#' monitor=c("m", "c", "precision"), n.chains=2,
#' method="rjags")
#' 
#' # Analyse standard plots of the results to assess convergence:
#' plot(results)
#' 
#' # Summary of the monitored variables:
#' results
#'
#' # For more details about possible methods see:
#' vignette('userguide', package='runjags')
#' }
#'

#' @param model either a relative or absolute path to a textfile (including the file extension) containing a model in the JAGS language and possibly monitored variable names, data and/or initial values, or a character string of the same.  No default.  See \code{\link{read.jagsfile}} for more details.

#' @param monitor a character vector of the names of variables to monitor.  No default.  The special node names 'deviance', 'pd', 'popt', 'dic', 'ped' and 'full.pd' are used to monitor the deviance, mean pD, mean pOpt, DIC, PED and full distribution of sum(pD) respectively.  Note that these monitored nodes (with the exception of 'deviance') require multiple chains within the same simulation, and won't appear as variables in the summary statistics or plots (but see \code{\link{extract}} for a way of extracting these from the returned object).

#' @param data a named list, data frame, environment, character string in the R dump format (see \code{\link{dump.format}}), or a function (with no arguments) returning one of these types.  If the model text contains inline #data# comments, then this argument specifies the list, data frame or environment in which to search first for these variables (the global environment is always searched last).  If the model text does not contain #data# comments, then the full list or data frame (but not environment) is included as data.  If the data argument is a character string, then any #data# comments in the model are ignored (with a warning). The default specifies the parent environment of the function call.

#' @param n.chains the number of chains to use with the simulation.  More chains will improve the sensitivity of the convergence diagnostic, but will cause the simulation to run more slowly (although this may be improved by using a method such as 'parallel', 'rjparallel' or 'snow').  The minimum (and default) number of chains is 2.

#' @param inits either a character vector with length equal to the number of chains the model will be run using, or a list of named lists representing names and corresponding values of inits for each chain, or  a function with either 1 argument representing the chain or no arguments.   If a vector, each element of the vector must be a character string in the  R dump format representing the initial values for that chain, or NA.  If not all initialising variables are specified, the unspecified variables are taken deterministically from the mean or mode of the prior distribution by JAGS.  Values left as NA result in all initial values for that chain being taken from the prior distribution.  The special variables '.RNG.seed', '.RNG.name', and '.RNG.state' are allowed for explicit control over random number generators in JAGS.  If a function is provided, the data is available inside the function as a named list 'data' - this may be useful for setting initial values that depend on the data.  Default NA.

#' @param runjags.object the model to be extended - the output of a run.jags (or autorun.jags or extend.jags etc) function, with class 'runjags'.  No default.

#' @param add.monitor a character vector of variables to add to the monitored variable list.  All previously monitored variables are automatically included - although see the 'drop.monitor' argument. Default no additional monitors.

#' @param drop.monitor a character vector of previously monitored variables to remove from the monitored variable list for the extended model. Default none.

#' @param drop.chain a numeric vector of chains to remove from the extended model. Default none.

#' @param combine a logical flag indicating if results from the new JAGS run should be combined with the previous chains.  Default TRUE if not adding or removing variables or chains, and FALSE otherwise.

#' @param startburnin the number of burnin iterations, NOT including the adaptive iterations to use for the initial pilot run of the chains.

#' @param startsample the total number of samples (including the chains supplied in runjags.object for autoextend.jags) on which to assess convergence, with a minimum of 4000.  If the runjags.object already contains this number of samples then convergence will be assessed on this object, otherwise the required number of additional samples will be obtained before combining the chains with the old chains. More samples will give a better chance of allowing the chain to converge, but will take longer to achieve.  Default 10000 iterations.

#' @param adapt the number of adaptive iterations to use at the start of each simulation.  For the rjags method this adaptation is only performed once and the model remains compiled, unless the repeatable.methods option is activated in \code{\link{runjags.options}}.  For all other methods adaptation is done every time the simulation is extended. Default 1000 iterations.

#' @param datalist deprecated argument.

#' @param initlist deprecated argument.

#' @param jags the system call or path for activating JAGS.  Default uses the option given in \code{\link{runjags.options}}.

#' @param silent.jags option to suppress output of the JAGS simulations.  Default uses the option given in \code{\link{runjags.options}}.

#' @param modules a character vector of external modules to be loaded into JAGS, either as the module name on its own or as the module name and status separated by a space, for example 'glm on'.

#' @param factories a character vector of factory modules to be loaded into JAGS.  Factories should be provided in the format '<facname> <factype> <status>' (where status is optional), for example: factories='mix::TemperedMix sampler on'.  You must also ensure that any required modules are also specified (in this case 'mix').

#' @param summarise should summary statistics be automatically calculated for the output chains?  Default TRUE (but see also ?runjags.options -> force.summary).

#' @param mutate either a function or a list with first element a function and remaining elements arguments to this function.  This can be used to add new variables to the posterior chains that are derived from the directly monitored variables in JAGS. This allows the variables to be summarised or extracted as part of the MCMC objects as if they had been calculated in JAGS, but without the computational or storage overheads associated with calculating them in JAGS directly.  The plot, summary and as.mcmc methods for runjags objects will automatically extract the mutated variables along with the directly monitored variables.  For an application to pairwise comparisons of different levels within fixed effects see \code{\link{contrasts.mcmc}}.

#' @param thin the thinning interval to be used in JAGS.  Increasing the thinning interval may reduce autocorrelation, and therefore reduce the number of samples required, but will increase the time required to run the simulation.  Using this option thinning is performed directly in JAGS, rather than on an existing MCMC object as with thin.sample. Default 1.

#' @param thin.sample option to thin the final MCMC chain(s) before calculating summary statistics and returning the chains.  Thinning very long chains reduces the size of the returned object.  If TRUE, the chain is thinned to as close to a minimum of startsample iterations as possible to ensure the chain length matches thin.sample.  A positive integer can also be specified as the desired chain length after thinning; the chains will be thinned to as close to this minimum value as possible. Default TRUE (thinned chains of length startsample returned).  This option does NOT carry out thinning in JAGS, therefore R must have enough available memory to hold the chains BEFORE thinning.  To avoid this problem use the 'thin' option instead.

#' @param raftery.options a named list which is passed as additional arguments to \code{\link[coda]{raftery.diag}}, or the logical FALSE to deactivate automatic run length calculation.  Default none (default arguments to raftery.diag are used).

#' @param crash.retry the number of times to re-attempt a simulation if the model returns an error.  Default 1 retry (simulation will be aborted after the second crash).

#' @param interactive option to allow the simulation to be interactive, in which case the user is asked if the simulation should be extended when run length and convergence calculations are performed and the extended simulation will take more than 1 minute.  The function will wait for a response before extending the simulations.  If FALSE, the simulation will be run until the chains have converged or until the next extension would extend the simulation beyond 'max.time'.  Default FALSE.

#' @param max.time the maximum time for which the function is allowed to extend the chains to improve convergence, as a character string including units or as an integer in which case units are taken as seconds.  Ignored if interactive=TRUE.  If the function thinks that the next simulation extension to improve convergence will result in a total time of greater than max.time, the extension is aborted.  The time per iteration is estimated from the first simulation.  Acceptable units include 'seconds', 'minutes', 'hours', 'days', 'weeks', or the first letter(s) of each.

#' @param tempdir option to use the temporary directory as specified by the system rather than creating files in the working directory.  Any files created in the temporary directory are removed when the function exits for any reason.  Default TRUE.

#' @param jags.refresh the refresh interval (in seconds) for monitoring JAGS output using the 'interactive' and 'parallel' methods (see the 'method' argument).  Longer refresh intervals will use slightly less processor time, but will make the simulation updates to be shown on the screen less frequently.  Reducing the refresh rate to every 10 or 30 seconds may be worthwhile for simulations taking several days to run.  Note that this will have no effect on the processor use of the simulations themselves.  Default 0.1 seconds.

#' @param batch.jags option to call JAGS in batch mode, rather than using input redirection.  On JAGS >= 3.0.0, this suppresses output of the status which may be useful in some situations.  Default TRUE if silent.jags is TRUE, or FALSE otherwise.

#' @param method the method with which to call JAGS; probably a character vector specifying one of 'rjags', 'simple', 'interruptible', 'parallel', 'rjparallel', or 'snow'. The 'rjags' and 'rjparallel' methods run JAGS using the rjags package, whereas other options do not require the rjags package and call JAGS as an external executable.  The advantage of the 'rjags' method is that the model will not need to be recompiled between successive calls to extend.jags, all other methods require a re-compilation (and adaptation if necessary) every time the model is extended.  The 'parallel' and 'interruptible' methods for Windows require XP Professional, Vista or later (or any Unix-alike).  For more information refer to the userguide vignette.

#' @param method.options a deprecated argument currently permitted for backwards compatibility, but this will be removed from a future version of runjags.  Pass these arguments directly to autorun.jags or autoextend.jags. 

#' @param ... summary parameters to be passed to \code{\link{add.summary}}, and/or additional options to control some methods including n.sims for parallel methods, cl for rjparallel and snow methods, remote.jags for snow, and by and progress.bar for the rjags method.
NULL

#' @rdname autorun.jags
autorun.jags <- function(model, monitor = NA, data=NA, n.chains=NA, inits = NA, startburnin = 4000, startsample = 10000, adapt=1000, datalist=NA, initlist=NA, jags = runjags.getOption('jagspath'), silent.jags = runjags.getOption('silent.jags'), modules=runjags.getOption('modules'), factories=runjags.getOption('factories'), summarise = TRUE, mutate = NA, thin = 1, thin.sample = FALSE, raftery.options = list(), crash.retry=1, interactive=FALSE, max.time=Inf, tempdir=runjags.getOption('tempdir'), jags.refresh=0.1, batch.jags=silent.jags, method=runjags.getOption('method'), method.options=list(), ...){
	
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
	
	res <- autoextend.jags(runjags.object=obj, add.monitor=character(0), drop.monitor=character(0), drop.chain=numeric(0), combine=FALSE, startburnin=startburnin, startsample=startsample, adapt=adapt, jags=jags, silent.jags=silent.jags, summarise=summarise, thin=thin, thin.sample=thin.sample, raftery.options=raftery.options, crash.retry=crash.retry, interactive=interactive, max.time=max.time, tempdir=tempdir, jags.refresh=jags.refresh, batch.jags=batch.jags, method=method, method.options=method.options, ...)
	
	return(res)
}
	

#' @rdname autorun.jags
autoextend.jags <- function(runjags.object, add.monitor=character(0), drop.monitor=character(0), drop.chain=numeric(0), combine=length(c(add.monitor,drop.monitor,drop.chain))==0, startburnin = 0, startsample = 10000, adapt=1000, jags = NA, silent.jags = NA, summarise = TRUE, thin = NA, thin.sample = FALSE, raftery.options = list(), crash.retry=1, interactive=FALSE, max.time=Inf, tempdir=runjags.getOption('tempdir'), jags.refresh=NA, batch.jags=NA, method=NA, method.options=NA, ...){

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
	
	# We may be passed some unevaluated function arguments so evaluate everything here:
	argnames <- names(formals(autoextend.jags))
	argnames <- argnames[argnames!='...']
	for(i in 1:length(argnames)){
		success <- try(assign(argnames[i], eval(get(argnames[i]))), silent=TRUE)		
		if(inherits(success, 'try-error')){
			stop(paste("object '", strsplit(as.character(success),split="'",fixed=TRUE)[[1]][2], "' not found", sep=""), call.=FALSE)
		}
	}
	
	# The summary parameters are checked by extend.jags and returned
	
	starttime <- Sys.time()
	
	# It's not clear which jags files to keep for autorun functions so disable:
	keep.jags.files <- FALSE
	
	method <- getrunjagsmethod(method)
	# Can't use background methods:
	if(method %in% c("background","bgparallel"))
		stop("The method specified to autorun.jags and autoextend.jags must run JAGS and wait for the results (ie the background method, and possibly other user specified methods, cannot be used)")

	# Check some autoextend specific stuff:
	if(length(drop.chain)>0){
		drop.chain <- unique(round(drop.chain))
		if(any(drop.chain>length(runjags.object$end.state)) | any(drop.chain < 1)) stop("Specified value(s) to drop.chain are invalid - please specify which chain(s) to drop by the chain number(s)")
		if(length(drop.chain)>(length(runjags.object$end.state)-2)) stop("Specified value(s) to drop.chains are invalid - a minimum of 2 chains must remain")
	}
	if(length(runjags.object$end.state) < 2)
		stop("The number of chains must be 2 or more so that convergence can be assessed",call.=FALSE)

	if(thin.sample==TRUE) thin.sample <- startsample
	if(thin.sample==FALSE) thin.sample <- Inf
	
	if(startsample!=0 && startsample < 4000)
		stop("A startsample of 4000 or more iterations (after thinning) is required to complete the Raftery and Lewis's diagnostic", call.=FALSE)
	
	# reftery.diag options are passed as a list or just FALSE (or list(FALSE)):
	if(identical(raftery.options, TRUE))
		raftery.options <- list()
	if(inherits(raftery.options,'list') && length(raftery.options)==1 && length(raftery.options[[1]])==1 && raftery.options[[1]]==FALSE){
		raftery.options <- FALSE
	}
	
	doraftery <- FALSE
	if(!identical(raftery.options, FALSE)){
		doraftery <- TRUE
		if(!is.list(raftery.options)) 
			stop("Options to raftery.diag must be provided as a named list")
		if(any(names(raftery.options)=="data")) 
			warning("The 'data' argument specified in raftery.options was ignored")
		raftery.args <- formals(raftery.diag)
		raftery.names <- names(raftery.args)
		if(length(raftery.options) > 0){
			if(is.null(names(raftery.options))) stop("Options to raftery.diag must be provided as a named list")
			raft.opt.names <- names(raftery.options)
			for(i in 1:length(raftery.options)){
				if(any(raft.opt.names[i]==raftery.names)){
					raftery.args[raft.opt.names[i]] <- raftery.options[i]
				}else{
					if(raft.opt.names[i]=="") stop("All arguments to raftery.diag must be named") else stop(paste(raft.opt.names[i], " is not a recognised argument to raftery.diag", sep=""))	
				}	
			}	

		}
		if(startsample!=0){	
			success <- try({
			raftery.args$data <- mcmc(1:startsample)
			class(raftery.args) <- "list"
			test <- do.call("raftery.diag", raftery.args)
			})	
			if(inherits(success, 'try-error')) stop("The arguments specified for raftery.diag are not valid")	
			if(test$resmatrix[1]=="Error") stop(paste("You need a startsample size of at least", test$resmatrix[2], "with the values of q, r and s specified for the raftery.options", sep=" "))
		}
	}
	
	
	# Get the maximum timeout:
	if(is.numeric(max.time)){
		max.time <- max.time #DEFAULT NOW SECONDS * 60
	}else{
		if(!is.character(max.time))
			stop("max.time must be either a numeric or character value")
		
		time.unit <- tolower(gsub('[^[:alpha:]]', '', max.time))
		if(time.unit=='secs')
			time.unit <- 'seconds'
		if(time.unit=='mins')
			time.unit <- 'minutes'		
		if(time.unit%in%c('hrs','hr'))
			time.unit <- 'hours'		
		possunits <- c('seconds','minutes','hours','days','weeks')
		matched.unit <- possunits[pmatch(time.unit, possunits)]
		if(is.na(matched.unit))
			stop(paste("Unrecognised unit of maximum time: '", time.unit, "'", sep=''))
		
		multiplier <- c(seconds=1, minutes=60, hours=60*60, days=24*60*60, weeks=24*60*60*7)[matched.unit]
		
		num.time <- suppressWarnings(as.numeric(gsub('[^[:digit:][:punct:]]', ' ', max.time)))
		if(is.na(num.time))
			stop(paste('Unable to extract a number from the value of "', max.time, '" specified to max.time', sep=''))
		
		max.time <- num.time * multiplier
		names(max.time) <- NULL
	}
	
	if(runjags.getOption('debug')){
		if(runjags.getOption('debug')>=10)
			print(max.time)
	}
	
	
	newlines <- if(silent.jags) "\n" else "\n\n"
	
	
	initialtimetaken <- runjags.object$timetaken

	# First run an extend.jags call with sample=0 - this deals with the add monitor and combine options, as well assummary.arg, method.options etc and checks JAGS etc etc - set silent.jags to TRUE though:
	runjags.object <- extend.jags(runjags.object, add.monitor=add.monitor, drop.monitor=drop.monitor, drop.chain=drop.chain, combine=combine && (runjags.object$sample>0), burnin = 0, sample = 0, adapt=0, noread.monitor = character(0), jags = jags, silent.jags = TRUE, summarise = FALSE, thin = thin, keep.jags.files = keep.jags.files, tempdir=tempdir, jags.refresh=jags.refresh, batch.jags=batch.jags, method=method, method.options=method.options, ...)
	# and then reset silentjags to what it should be:
	runjags.object$method.options$silent.jags <- silent.jags

	# Catch where we want a model with no updates and return:
	if(startsample==0)
		return(runjags.object)
	
	summaryargs <- runjags.object$summary.pars
	psrf.target <- summaryargs$psrf.target
	method.options <- runjags.object$method.options
	
	swcat("\nAuto-run JAGS",newlines,sep="")	

	# This function call gave the necessary initial value and data etc warnings, so make sure we don't get them again:
	currentinitwarn <- runjags.getOption('inits.warning')
	currentdatawarn <- runjags.getOption('nodata.warning')
	currentrecover <- runjags.getOption('partial.import')
	on.exit(runjags.options(inits.warning=currentinitwarn, nodata.warning=currentdatawarn, partial.import=currentrecover))
	runjags.options(inits.warning=FALSE)
	runjags.options(nodata.warning=FALSE)
	runjags.options(partial.import=FALSE)
	
	if(startsample > runjags.object$sample){
		
		swcat("Running a pilot chain...\n")
		
		# runjags.object$sample and $burnin may be 0:
		initialsample <- startsample - runjags.object$sample
		
		# Run for some more iterations:
		additional <- extend.jags(runjags.object, combine=runjags.object$sample>0, burnin=startburnin, sample=initialsample, adapt=adapt, jags = jags, silent.jags = silent.jags, summarise = FALSE, thin = thin, keep.jags.files = keep.jags.files, tempdir=tempdir, jags.refresh=jags.refresh, batch.jags=batch.jags)
		if(class(additional)=='runjagsbginfo') stop("The method specified to autorun.jags and autoextend.jags must run JAGS and wait for the results (ie the background method, and possibly other user specified methods, cannot be used)")
				
		if(niter(additional$mcmc) < initialsample){
			repeat{
				time.taken <- timestring(starttime, Sys.time(), units="secs", show.units=FALSE)
				if(time.taken > max.time | crash.retry==0){
					stop("The simulation exceeded the number of crashes allowed by crash.retry and so was aborted", call.=FALSE)
				}
				swcat("\nThe simulation crashed; retrying...",newlines,sep="")			
				crash.retry <- crash.retry - 1

				additional <- extend.jags(runjags.object, combine=TRUE, burnin=startburnin, sample=initialsample, adapt=adapt, jags = jags, silent.jags = silent.jags, summarise = FALSE, thin = thin, keep.jags.files = keep.jags.files, tempdir=tempdir, jags.refresh=jags.refresh, batch.jags=batch.jags)

				if(niter(additional$mcmc) == initialsample) break
			}
		}
		
		swcat("\n")
	
	}else{
		additional <- runjags.object
	}
	firsttimetaken = time.taken <- timestring(starttime, Sys.time(), units="secs", show.units=FALSE)

	swcat("Calculating the Gelman-Rubin statistic for ", nvar(additional$mcmc), " variables....\n", sep="")
	
	suppressWarnings(success <- try(convergence <- safe.gelman.diag(normalise.mcmcfun(additional$mcmc, normalise=summaryargs$normalise.mcmc, warn=FALSE, remove.nonstochastic = TRUE)$mcmc, transform=FALSE, autoburnin=FALSE), silent=TRUE))
	if(inherits(success, 'try-error')){
		stop("An error occured while calculating the Gelman-Rubin statistic.  Check that different chains have not been given the same starting values and random seeds, and that there is at least one monitored stochastic variable.",call.=FALSE)
	}
	
	convergence <- c(convergence, psrf.target=psrf.target)
	class(convergence) <- "gelmanwithtarget"
	
	n.params <- nrow(convergence$psrf)
	n.iters <- niter(additional$mcmc)
	
	if(n.params==1) convergence$mpsrf <- convergence$psrf[1,1]
	unconverged <- 0

	for(j in 1:n.params){
		param.conv <- convergence$psrf[j, 1]
		if(!is.na(param.conv)){
			if(param.conv > psrf.target){
				unconverged <- unconverged + 1
			}
		}else{
			warning(paste("The Gelman-Rubin statistic for '", varnames(additional$mcmc)[j], "' could not be calculated", sep=""))
		}
	}

	if(unconverged > 0){
		if(!is.numeric(convergence$mpsrf)){
			mpsrfstring <- " (Unable to calculate the multi-variate psrf)"
		}else{
			mpsrfstring <- paste(" (multi-variate psrf = ", round(convergence$mpsrf, digits=3), ")", sep="")
		}
			
		swcat("The Gelman-Rubin statistic was above ", psrf.target, " for ", unconverged, " parameter", if(unconverged>1) "s", " after ", additional$burnin+additional$sample, " iterations taking ", timestring(additional$timetaken), " ", mpsrfstring, ".  This may indicate poor convergence.\n", sep="")
		time.taken <- timestring(starttime, Sys.time(), units="secs", show.units=FALSE)
		stop <- time.taken > max.time
			
		if(interactive){
			time.taken <- timestring(starttime, Sys.time(), units="secs", show.units=FALSE)
			stop <- !ask("Extend the simulation to attempt to improve convergence?")
			starttime <- Sys.time()-time.taken
		}
		if(stop){
			
			if(!interactive)
				swcat("Maximum time limit exceeded; chains unconverged.  Try restarting the simulation using the end state values of the chains provided.\n")
				
			swcat("Calculating autocorrelation and summary statistics...\n")
			
			timetaken <- (difftime(Sys.time(), starttime, units='secs') + initialtimetaken)
	
			if(niter(additional$mcmc)>thin.sample){
				additional$mcmc <- combine.mcmc(additional$mcmc, collapse.chains=FALSE, return.samples=thin.sample)
				if(!identical(additional$pd, NA)) additional$pd <- combine.mcmc(additional$pd, collapse.chains=FALSE, return.samples=thin.sample)
			}
	
			combinedoutput <- makerunjagsobject(additional, summarise=summarise, summaryargs=summaryargs, burnin=additional$burnin, sample=niter(additional$mcmc), thin=thin, model=runjags.object$model, data=runjags.object$data, monitor=runjags.object$monitor, noread.monitor=runjags.object$noread.monitor, modules=runjags.object$modules, factories=runjags.object$factories, response=runjags.object$response, residual=runjags.object$residual, fitted=runjags.object$fitted, method=runjags.object$method, method.options=runjags.object$method.options, timetaken=timetaken)
				
			swcat("Returning UNCONVERGED simulation results\n\n")
			return(combinedoutput)
		}
		
		finishconv <- FALSE
		swcat("Extending the simulation to attempt to improve convergence...\n")
						
		while(!finishconv){
			
			# Run for some more iterations:
			extended <- extend.jags(additional, combine=FALSE, burnin=0, sample=startsample, adapt=adapt, jags = jags, silent.jags = silent.jags, summarise = FALSE, thin = thin, keep.jags.files = keep.jags.files, tempdir=tempdir, jags.refresh=jags.refresh, batch.jags=batch.jags)
					
			if(niter(extended$mcmc) < startsample){
				repeat{
					time.taken <- timestring(starttime, Sys.time(), units="secs", show.units=FALSE)
					if(time.taken > max.time | crash.retry==0){
						stop("The simulation exceeded the number of crashes allowed by crash.retry and so was aborted", call.=FALSE)
					}
					swcat("\nThe simulation crashed; retrying...",newlines,sep="")			
					crash.retry <- crash.retry - 1
					
					extended <- extend.jags(additional, combine=FALSE, burnin=0, sample=startsample, adapt=adapt, jags = jags, silent.jags = silent.jags, summarise = FALSE, thin = thin, keep.jags.files = keep.jags.files, tempdir=tempdir, jags.refresh=jags.refresh, batch.jags=batch.jags)
					
					if(niter(extended$mcmc) == startsample) break
				}
			}
			
			additional <- extended

			swcat("Calculating the Gelman-Rubin statistic for ", nvar(additional$mcmc), " variables....\n", sep="")
			suppressWarnings(success <- try(convergence <- safe.gelman.diag(normalise.mcmcfun(additional$mcmc, normalise=summaryargs$normalise.mcmc, warn=FALSE, remove.nonstochastic = TRUE)$mcmc, transform=FALSE, autoburnin=FALSE), silent=TRUE))
			if(inherits(success, 'try-error')){
				stop("An error occured while calculating the Gelman-Rubin statistic.  Check that different chains have not been given the same starting values and random seeds, and that there is at least one stochastic monitored variable.", call.=FALSE)
			}
				
			convergence <- c(convergence, psrf.target=psrf.target)
			class(convergence) <- "gelmanwithtarget"
				
			n.params <- nrow(convergence$psrf)				
				
			unconverged <- 0
			if(n.params==1) convergence$mpsrf <- convergence$psrf[1,1]
			
			for(j in 1:n.params){
				param.conv <- convergence$psrf[j, 1]
				if(!is.na(param.conv)){
					if(param.conv > psrf.target){
						unconverged <- unconverged + 1
					}
				}else{
					warning(paste("The Gelman-Rubin statistic for '", varnames(additional$mcmc)[j], "' could not be calculated", sep=""))
				}	
			}
									
			if(is.numeric(convergence$mpsrf)){
				mpsrfstring <- " (Unable to calculate the multi-variate psrf)"
			}else{
				mpsrfstring <- paste(" (multi-variate psrf = ", round(convergence$mpsrf, digits=3), ")", sep="")
			}
				
			if(unconverged > 0){
				swcat("The Gelman-Rubin statistic was still above ", psrf.target, " for ", unconverged, " parameter", if(unconverged>1) "s", " after ", additional$burnin+additional$sample, " iterations taking ", timestring(additional$timetaken), " ", mpsrfstring, ".\n", sep="")
				
				stop <- difftime(Sys.time(), starttime, units='secs') > max.time
				
				if(interactive){
					stop <- !ask("Extend the simulation to attempt to improve convergence?")
					starttime <- Sys.time()
				}
				if(stop){
					
					if(!interactive)
						swcat("Maximum time limit exceeded; chains still unconverged.  Try restarting the simulation using the end state values of the chains provided.\n")
					
					swcat("Calculating autocorrelation and summary statistics...\n")
		
					timetaken <- (difftime(Sys.time(), starttime, units='secs') + initialtimetaken)

					if(niter(additional$mcmc)>thin.sample){
						additional$mcmc <- combine.mcmc(additional$mcmc, collapse.chains=FALSE, return.samples=thin.sample)
						if(!identical(additional$pd, NA)) additional$pd <- combine.mcmc(additional$pd, collapse.chains=FALSE, return.samples=thin.sample)
					}

					combinedoutput <- makerunjagsobject(additional, summarise=summarise, summaryargs=summaryargs, burnin=additional$burnin, sample=niter(additional$mcmc), thin=thin, model=runjags.object$model, data=runjags.object$data, monitor=runjags.object$monitor, noread.monitor=runjags.object$noread.monitor, modules=runjags.object$modules, factories=runjags.object$factories, response=runjags.object$response, residual=runjags.object$residual, fitted=runjags.object$fitted, method=runjags.object$method, method.options=runjags.object$method.options, timetaken=timetaken)
			
					swcat("Returning UNCONVERGED simulation results\n\n")
					return(combinedoutput)
												
				}else{
					swcat("Extending the simulation to attempt to improve convergence...\n")
				}
					
			}else{
				swcat(paste("The Gelman-Rubin statistic is now below ", psrf.target, " for all parameters\n", sep=""))
				finishconv <- TRUE
			}
		}  # End of the while loop
	
	}else{
		swcat(paste("The Gelman-Rubin statistic is below ", psrf.target, " for all parameters\n", sep=""))
	}	
	
	
	moreupdates <- 0
	if(doraftery){
		swcat("\nCalculating the necessary sample length based on the Raftery and Lewis's diagnostic...\n")
			
		success <- try({
		raftery.args$data <- normalise.mcmcfun(additional$mcmc, normalise=FALSE, warn=FALSE, remove.nonstochastic = TRUE)$mcmc
		class(raftery.args) <- "list"
		raftery <- do.call("raftery.diag", raftery.args)
		})
	
		if(inherits(success, 'try-error')) stop("An error occured while calculating the Raftery and Lewis's diagnostic",call.=FALSE)
		if(raftery[[1]]$resmatrix[1]=="error") stop("Error", "An error occured while calculating the Raftery and Lewis diagnostic",call.=FALSE)
	
		# to correct for monitoring arrays and non-stochastic nodes:
		newmonitor <- dimnames(raftery[[1]]$resmatrix)[[1]]
		n.chains <- length(additional$end.state)
	
		dependance = burnin = sample <- matrix(ncol=n.chains, nrow=length(newmonitor), dimnames=list(dimnames(raftery[[1]]$resmatrix)[[1]], 1:n.chains))
	
		for(i in 1:n.chains){	
			dependance[,i] <- raftery[[i]]$resmatrix[,"I"]
			burnin[,i] <- raftery[[i]]$resmatrix[,"M"]
			sample[,i] <- raftery[[i]]$resmatrix[,"N"]
		}
	
		dependancethreshold <- 3
	
	#	if(any(dependance > dependancethreshold) & killautocorr==FALSE){
	#		swcat("IMPORTANT:  The sample size of monitored node(s) '", paste(dimnames(dependance)[[1]][apply(dependance, 1, function(x) if(any(x>dependancethreshold)) return(TRUE) else return(FALSE))], collapse="' & '"), "' have a high autocorrelation dependance in chain(s) ", paste(seq(1, n.chains)[apply(dependance, 2, function(x) if(any(x>dependancethreshold)) return(TRUE) else return(FALSE))], collapse= " & "), ".  Re-running the model with a different formulation or better initial values may help to reduce autocorrelation.\n", sep="")
	#	}
	
		#### raftery.diag takes account of the chain thinning already, so we need to multiply the number of iterations done by the thinning to work out what we have left:
		totalunthinnedperchain <- max(sample)/n.chains
		moreupdates <- max(totalunthinnedperchain - (thin*niter(additional$mcmc)), 0) / thin
		moreupdates <- ceiling(moreupdates)
	
		if(runjags.getOption('debug')){
			if(runjags.getOption('debug')>=10)
				print(raftery)
		
			swcat('Raftery diag:  maxsample=', max(sample), '/', thin, '=', max(sample)/thin, ', current=', niter(additional$mcmc), '*', n.chains, '=', niter(additional$mcmc)*n.chains, ', more=', moreupdates, '*', n.chains, '\n', sep='')

		}
	}
	
	if(moreupdates > 0){
		
		# There must be a better way to do this...
		success <- try({
			testmatrix <- matrix(nrow=(max(sample)/thin)*n.chains, ncol=length(varnames(additional$mcmc)))
			rm(testmatrix)
		})
		if(inherits(success, 'try-error')){
			stop(paste("The model needs to be run for a further ", moreupdates, " iterations.  This would create a vector too large to be read into R.  Try re-parameterising the model to reduce autocorrelation, using the thin option to reduce autocorrelation, or monitoring less variables.  The simulation can be extended using the return value provided.", sep=""), call.=FALSE)
		}
		
		swcat("The model will need to be run for a further ", moreupdates, " updates.  This will take approximately ", timestring((firsttimetaken*moreupdates/(startsample+startburnin))), ".\n", sep="")
		if(interactive & (timestring((time.taken*moreupdates/(startsample+startburnin)), units="s", show.units=FALSE)>60)){
			time.taken <- timestring(starttime, Sys.time(), units="secs", show.units=FALSE)
			if(!ask("Continue with the simulation?")){
					
				swcat("Calculating autocorrelation and summary statistics...\n")
	
				timetaken <- (difftime(Sys.time(), starttime, units='secs') + initialtimetaken)

				if(niter(additional$mcmc)>thin.sample){
					additional$mcmc <- combine.mcmc(additional$mcmc, collapse.chains=FALSE, return.samples=thin.sample)
					if(!identical(additional$pd, NA)) additional$pd <- combine.mcmc(additional$pd, collapse.chains=FALSE, return.samples=thin.sample)
				}

				combinedoutput <- makerunjagsobject(additional, summarise=summarise, summaryargs=summaryargs, burnin=additional$burnin, sample=niter(additional$mcmc), thin=thin, model=runjags.object$model, data=runjags.object$data, monitor=runjags.object$monitor, noread.monitor=runjags.object$noread.monitor, modules=runjags.object$modules, factories=runjags.object$factories, response=runjags.object$response, residual=runjags.object$residual, fitted=runjags.object$fitted, method=runjags.object$method, method.options=runjags.object$method.options, timetaken=timetaken)
		
				swcat("Simulation aborted\n\n")
				return(combinedoutput)
					
			}
		}
		
		swcat("\n")
		
		# Run for some more iterations:
		extended <- extend.jags(additional, combine=TRUE, burnin=0, sample=moreupdates, adapt=adapt, jags = jags, silent.jags = silent.jags, summarise = FALSE, thin = thin, keep.jags.files = keep.jags.files, tempdir=tempdir, jags.refresh=jags.refresh, batch.jags=batch.jags)
		
		if(niter(extended$mcmc) < moreupdates){
			repeat{
				time.taken <- timestring(starttime, Sys.time(), units="secs", show.units=FALSE)
				if(time.taken > max.time | crash.retry==0){
					stop("The simulation exceeded the number of crashes allowed by crash.retry and so was aborted", call.=FALSE)
				}
				swcat("\nThe simulation crashed; retrying...",newlines,sep="")			
				crash.retry <- crash.retry - 1
				
				extended <- extend.jags(additional, combine=TRUE, burnin=0, sample=moreupdates, adapt=adapt, jags = jags, silent.jags = silent.jags, summarise = FALSE, thin = thin, keep.jags.files = keep.jags.files, tempdir=tempdir, jags.refresh=jags.refresh, batch.jags=batch.jags)
								
				if(niter(extended$mcmc) == moreupdates) break
			}
		}
		
		additional <- extended
	}
	if(doraftery)
		swcat("Indicated sample length achieved\n")
	
	# Account for thin.sample:
	if(niter(additional$mcmc)>thin.sample){
		additional$mcmc <- combine.mcmc(additional$mcmc, collapse.chains=FALSE, return.samples=thin.sample)
		if(!identical(additional$pd, NA)) additional$pd <- combine.mcmc(additional$pd, collapse.chains=FALSE, return.samples=thin.sample)
	}
	
	timetaken <- (difftime(Sys.time(), starttime) + initialtimetaken)

	combinedoutput <- makerunjagsobject(additional, summarise=summarise, summaryargs=summaryargs, burnin=additional$burnin, sample=niter(additional$mcmc), thin=thin, model=runjags.object$model, data=runjags.object$data, monitor=runjags.object$monitor, noread.monitor=runjags.object$noread.monitor, modules=runjags.object$modules, factories=runjags.object$factories, response=runjags.object$response, residual=runjags.object$residual, fitted=runjags.object$fitted, method=runjags.object$method, method.options=runjags.object$method.options, timetaken=timetaken)

	swcat("Auto-run JAGS complete\n\n")
	
	return(combinedoutput)
	
}

autorun.JAGS <- autorun.jags
autoextend.JAGS <- autoextend.jags
