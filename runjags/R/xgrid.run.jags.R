xgrid.run.jags <- function(model, max.threads=Inf, JAGSversion=">=2.0.0", email=NA, profiling=TRUE, cpuarch=NA, minosversion=NA, queueforserver=FALSE, hostnode=NA, forcehost=FALSE, ramrequired=10, jobname=NA, cleanup=TRUE, showprofiles=FALSE, jagspath='/usr/local/bin/jags', mgridpath=system.file("xgrid", "mgrid.sh", package="runjags"), hostname=Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), password=Sys.getenv("XGRID_CONTROLLER_PASSWORD"), ...){
	
	if(max.threads < 1) stop("The maximum number of threads must be greater than 0")
		
	passed <- c(list(...), model=model)
	if(any(names(passed)=="method") || any(names(passed)=="method.options" )) stop("Cannot specify 'method' or 'method.options' arguments to xgrid functions")
	
	# Translate monitor.pd/popt/pd.i/deviance (legacy code):
	monitor.deviance <- any(names(passed)=="monitor.deviance") && passed$monitor.deviance
	monitor.pd <- any(names(passed)=="monitor.pd") && passed$monitor.pd
	monitor.popt <- any(names(passed)=="monitor.popt") && passed$monitor.popt
	monitor.pd.i <- any(names(passed)=="monitor.pd.i") && passed$monitor.pd.i
	if(any(names(passed)=="check.conv")){
		warning("Use of the 'check.conv' argument is deprecated - use 'summarise' to achieve some of the same function", call.=FALSE)
		if(!any(names(passed)=="summarise")) passed <- c(passed, list(summarise=passed$check.conv))
	}

	if(any(names(passed)=="monitor")){
		passed$monitor <- c(passed$monitor, if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")
	}else{
		if(any(c(monitor.deviance,monitor.pd,monitor.popt,monitor.pd.i))){
			passed <- c(passed, list(c(if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")))
			warning("Use of the 'monitor.deviance', 'monitor.pd', 'monitor.popt' and 'monitor.pdi' arguments is deprecated - add the special variables 'dic', 'deviance', 'pd', 'popt' or 'pd.i' to the 'monitor' argument instead", call.=FALSE)
		}
	}
	
	# Get default argument values from specified functions for anything not given to this function:
 	arglist <- getargs(c('run.jags','add.summary'), passed)	
	
 	args <- arglist[which(names(arglist) %in% names(formals(setup.jagsfile)))]
	args$method <- if(any(.packages()=="rjags")) "rjags" else "simple"
	if(!any(names(args)=='data'))
		args$data <- parent.frame()
	if(!any(names(args)=='inits'))
		args$inits <- parent.frame()
 	obj <- do.call("setup.jagsfile", args=args)

	arglist <- arglist[names(arglist)!="method" & names(arglist)!="method.options"]
	
	# Get xgrid method options:
	method.options <- setup.xgrid(separate=max.threads>1,JAGSversion=JAGSversion, Rversion="", packages=list(), artfun=function() writeLines("1"), email=email, profiling=profiling, cpuarch=cpuarch, minosversion=minosversion, queueforserver=queueforserver, hostnode=hostnode, forcehost=forcehost, ramrequired=ramrequired, jobname=jobname, cleanup=cleanup, showprofiles=showprofiles, jagspath=jagspath, Rpath='/usr/bin/R', Rbuild='64', max.filesize="1GB", mgridpath=mgridpath, hostname=hostname, password=password, submitandstop=FALSE, jagsrun=TRUE)
	method.options$max.threads <- max.threads
	method.options <- method.options[c("jobname", "command", "customart", "jagspath", "submitandstop", "max.threads", "mgridpath", "hostname", "password")]
	
	if(max.threads>1 && any(obj$monitor %in% c("dic", "pd", "popt", "pd.i"))){
		warning("Setting max.threads to 1 to allow calculation of dic/pd/popt/pd.i")
		max.threads <- 1
	}
	method.options$max.threads <- max.threads
	
	# Slight hack - we're not calling actually run.jags so add combine=FALSE to the argument list:
	arglist <- c(arglist, list(runjags.object=obj, combine=FALSE, method="xgrid", method.options=method.options)) 
	
	if(!any(names(arglist)=="burnin")) arglist$burnin <- 5000
	
 	args <- arglist[which(names(arglist) %in% c(names(formals(extend.jags)), names(formals(add.summary))))]
	
 	res <- do.call("extend.jags", args=args)
	
	warning("As of runjags version 2.0, the xgrid.xxx.jags functions are untested. If you are still using them, please email the package maintainer - otherwise they are likely to be removed from a future release.")
	return(res)
}


xgrid.autorun.jags <- function(model, max.threads=Inf, JAGSversion=">=2.0.0", email=NA, profiling=TRUE, cpuarch=NA, minosversion=NA, queueforserver=FALSE, hostnode=NA, forcehost=FALSE, ramrequired=10, jobname=NA, cleanup=TRUE, showprofiles=FALSE, jagspath='/usr/local/bin/jags', mgridpath=system.file("xgrid", "mgrid.sh", package="runjags"), hostname=Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), password=Sys.getenv("XGRID_CONTROLLER_PASSWORD"), ...){
	
	if(max.threads < 1) stop("The maximum number of threads must be greater than 0")
		
	passed <- c(list(...), model=model)
	if(any(names(passed)=="method") || any(names(passed)=="method.options" )) stop("Cannot specify 'method' or 'method.options' arguments to xgrid functions")
	
	# Translate monitor.pd/popt/pd.i/deviance (legacy code):
	monitor.deviance <- any(names(passed)=="monitor.deviance") && passed$monitor.deviance
	monitor.pd <- any(names(passed)=="monitor.pd") && passed$monitor.pd
	monitor.popt <- any(names(passed)=="monitor.popt") && passed$monitor.popt
	monitor.pd.i <- any(names(passed)=="monitor.pd.i") && passed$monitor.pd.i
	
	if(any(names(passed)=="monitor")){
		passed$monitor <- c(passed$monitor, if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")
	}else{
		if(any(c(monitor.deviance,monitor.pd,monitor.popt,monitor.pd.i))){
			passed <- c(passed, list(c(if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")))
			warning("Use of the 'monitor.deviance', 'monitor.pd', 'monitor.popt' and 'monitor.pdi' arguments is deprecated - add the special variables 'dic', 'deviance', 'pd', 'popt' or 'pd.i' to the 'monitor' argument instead", call.=FALSE)
		}
	}
	
	# Get default argument values from specified functions for anything not given to this function:
 	arglist <- getargs(c('autorun.jags','add.summary'), passed)	
	
 	args <- arglist[which(names(arglist) %in% names(formals(setup.jagsfile)))]
	args$method <- if(any(.packages()=="rjags")) "rjags" else "simple"
	if(!any(names(args)=='data'))
		args$data <- parent.frame()
	if(!any(names(args)=='inits'))
		args$inits <- parent.frame()
	
 	obj <- do.call("setup.jagsfile", args=args)
	
	arglist <- arglist[names(arglist)!="method" & names(arglist)!="method.options"]
	
	# Get xgrid method options:
	method.options <- setup.xgrid(separate=max.threads>1,JAGSversion=JAGSversion, Rversion="", packages=list(), artfun=function() writeLines("1"), email=email, profiling=profiling, cpuarch=cpuarch, minosversion=minosversion, queueforserver=queueforserver, hostnode=hostnode, forcehost=forcehost, ramrequired=ramrequired, jobname=jobname, cleanup=cleanup, showprofiles=showprofiles, jagspath=jagspath, Rpath='/usr/bin/R', Rbuild='64', max.filesize="1GB", mgridpath=mgridpath, hostname=hostname, password=password, submitandstop=FALSE, jagsrun=TRUE)
	method.options$max.threads <- max.threads
	method.options <- method.options[c("jobname", "command", "customart", "jagspath", "submitandstop", "max.threads", "mgridpath", "hostname", "password")]
	
	if(max.threads>1 && any(obj$monitor %in% c("dic", "pd", "popt", "pd.i"))){
		warning("Setting max.threads to 1 to allow calculation of dic/pd/popt/pd.i")
		max.threads <- 1
	}
	method.options$max.threads <- max.threads
	
	# Slight hack - we're not calling actually run.jags so add combine=FALSE to the argument list:
	arglist <- c(arglist, list(runjags.object=obj, combine=FALSE, method="xgrid", method.options=method.options))
	if(!any(names(arglist)=="startburnin")) arglist$startburnin <- 5000
	
 	args <- arglist[which(names(arglist) %in% c(names(formals(autoextend.jags)), names(formals(add.summary))))]
	
 	res <- do.call("autoextend.jags", args=args)
	
	warning("As of runjags version 2.0, the xgrid.xxx.jags functions are untested. If you are still using them, please email the package maintainer - otherwise they are likely to be removed from a future release.")
	return(res)
}

xgrid.extend.jags <- function(runjags.object, max.threads=Inf, JAGSversion=">=2.0.0", email=NA, profiling=TRUE, cpuarch=NA, minosversion=NA, queueforserver=FALSE, hostnode=NA, forcehost=FALSE, ramrequired=10, jobname=NA, cleanup=TRUE, showprofiles=FALSE, jagspath='/usr/local/bin/jags', mgridpath=system.file("xgrid", "mgrid.sh", package="runjags"), hostname=Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), password=Sys.getenv("XGRID_CONTROLLER_PASSWORD"), ...){
	
	checkvalidrunjagsobject(runjags.object)

	if(max.threads < 1) stop("The maximum number of threads must be greater than 0")
		
	# We are only passing straight through to extend.jags so don't really need all this but may as well  keep it in:
	
	passed <- c(list(...), list(x=runjags.object))
	if(any(names(passed)=="method") || any(names(passed)=="method.options" )) stop("Cannot specify 'method' or 'method.options' arguments to xgrid functions")
	
	# Translate monitor.pd/popt/pd.i/deviance (legacy code):
	monitor.deviance <- any(names(passed)=="monitor.deviance") && passed$monitor.deviance
	monitor.pd <- any(names(passed)=="monitor.pd") && passed$monitor.pd
	monitor.popt <- any(names(passed)=="monitor.popt") && passed$monitor.popt
	monitor.pd.i <- any(names(passed)=="monitor.pd.i") && passed$monitor.pd.i
	if(any(names(passed)=="check.conv")){
		warning("Use of the 'check.conv' argument is deprecated - use 'summarise' to achieve some of the same function", call.=FALSE)
		if(!any(names(passed)=="summarise")) passed <- c(passed, list(summarise=passed$check.conv))
	}
	
	if(any(names(passed)=="monitor")){
		passed$monitor <- c(passed$monitor, if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")
	}else{
		if(any(c(monitor.deviance,monitor.pd,monitor.popt,monitor.pd.i))){
			passed <- c(passed, list(c(if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")))
			warning("Use of the 'monitor.deviance', 'monitor.pd', 'monitor.popt' and 'monitor.pdi' arguments is deprecated - add the special variables 'dic', 'deviance', 'pd', 'popt' or 'pd.i' to the 'monitor' argument instead", call.=FALSE)
		}
	}
	
	# Get default argument values from specified functions for anything not given to this function:
 	arglist <- getargs(c('extend.jags','add.summary'), passed)	
 	arglist <- arglist[names(arglist)!="method" & names(arglist)!="method.options"]
	
	# Get xgrid method options:
	method.options <- setup.xgrid(separate=max.threads>1,JAGSversion=JAGSversion, Rversion="", packages=list(), artfun=function() writeLines("1"), email=email, profiling=profiling, cpuarch=cpuarch, minosversion=minosversion, queueforserver=queueforserver, hostnode=hostnode, forcehost=forcehost, ramrequired=ramrequired, jobname=jobname, cleanup=cleanup, showprofiles=showprofiles, jagspath=jagspath, Rpath='/usr/bin/R', Rbuild='64', max.filesize="1GB", mgridpath=mgridpath, hostname=hostname, password=password, submitandstop=FALSE, jagsrun=TRUE)
	method.options$max.threads <- max.threads
	method.options <- method.options[c("jobname", "command", "customart", "jagspath", "submitandstop", "max.threads", "mgridpath", "hostname", "password")]
	
	if(max.threads>1 && any(runjags.object$monitor %in% c("dic", "pd", "popt", "pd.i"))){
		warning("Setting max.threads to 1 to allow calculation of dic/pd/popt/pd.i")
		max.threads <- 1
	}
	method.options$max.threads <- max.threads
	
	arglist <- c(arglist, list(method="xgrid", method.options=method.options))
	args <- arglist[which(names(arglist) %in% c(names(formals(extend.jags)), names(formals(add.summary))))]
	
 	res <- do.call("extend.jags", args=args)
	
	warning("As of runjags version 2.0, the xgrid.xxx.jags functions are untested. If you are still using them, please email the package maintainer - otherwise they are likely to be removed from a future release.")
	return(res)
}

xgrid.autoextend.jags <- function(runjags.object, max.threads=Inf, JAGSversion=">=2.0.0", email=NA, profiling=TRUE, cpuarch=NA, minosversion=NA, queueforserver=FALSE, hostnode=NA, forcehost=FALSE, ramrequired=10, jobname=NA, cleanup=TRUE, showprofiles=FALSE, jagspath='/usr/local/bin/jags', mgridpath=system.file("xgrid", "mgrid.sh", package="runjags"), hostname=Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), password=Sys.getenv("XGRID_CONTROLLER_PASSWORD"), ...){
	
	checkvalidrunjagsobject(runjags.object)

	if(max.threads < 1) stop("The maximum number of threads must be greater than 0")
		
	# We are only passing straight through to autoextend.jags so don't really need all this but may as well keep it in:

	passed <- c(list(...), list(x=runjags.object))
	if(any(names(passed)=="method") || any(names(passed)=="method.options" )) stop("Cannot specify 'method' or 'method.options' arguments to xgrid functions")
	
	# Translate monitor.pd/popt/pd.i/deviance (legacy code):
	monitor.deviance <- any(names(passed)=="monitor.deviance") && passed$monitor.deviance
	monitor.pd <- any(names(passed)=="monitor.pd") && passed$monitor.pd
	monitor.popt <- any(names(passed)=="monitor.popt") && passed$monitor.popt
	monitor.pd.i <- any(names(passed)=="monitor.pd.i") && passed$monitor.pd.i
	
	if(any(names(passed)=="monitor")){
		passed$monitor <- c(passed$monitor, if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")
	}else{
		if(any(c(monitor.deviance,monitor.pd,monitor.popt,monitor.pd.i))){
			passed <- c(passed, list(c(if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")))
			warning("Use of the 'monitor.deviance', 'monitor.pd', 'monitor.popt' and 'monitor.pdi' arguments is deprecated - add the special variables 'dic', 'deviance', 'pd', 'popt' or 'pd.i' to the 'monitor' argument instead", call.=FALSE)
		}
	}
	
	# Get default argument values from specified functions for anything not given to this function:
 	arglist <- getargs(c('autoextend.jags','add.summary'), passed)	
	
	arglist <- arglist[names(arglist)!="method" & names(arglist)!="method.options"]
	
	# Get xgrid method options:
	method.options <- setup.xgrid(separate=max.threads>1,JAGSversion=JAGSversion, Rversion="", packages=list(), artfun=function() writeLines("1"), email=email, profiling=profiling, cpuarch=cpuarch, minosversion=minosversion, queueforserver=queueforserver, hostnode=hostnode, forcehost=forcehost, ramrequired=ramrequired, jobname=jobname, cleanup=cleanup, showprofiles=showprofiles, jagspath=jagspath, Rpath='/usr/bin/R', Rbuild='64', max.filesize="1GB", mgridpath=mgridpath, hostname=hostname, password=password, submitandstop=FALSE, jagsrun=TRUE)
	method.options$max.threads <- max.threads
	method.options <- method.options[c("jobname", "command", "customart", "jagspath", "submitandstop", "max.threads", "mgridpath", "hostname", "password")]
	
	if(max.threads>1 && any(runjags.object$monitor %in% c("dic", "pd", "popt", "pd.i"))){
		warning("Setting max.threads to 1 to allow calculation of dic/pd/popt/pd.i")
		max.threads <- 1
	}
	method.options$max.threads <- max.threads
	
	arglist <- c(arglist, list(method="xgrid", method.options=method.options))
 	args <- arglist[which(names(arglist) %in% c(names(formals(autoextend.jags)), names(formals(add.summary))))]
	
 	res <- do.call("autoextend.jags", args=args)
	
	warning("As of runjags version 2.0, the xgrid.xxx.jags functions are untested. If you are still using them, please email the package maintainer - otherwise they are likely to be removed from a future release.")
	return(res)
}

xgrid.submit.jags <- function(model, max.threads=Inf, JAGSversion=">=2.0.0", email=NA, profiling=TRUE, cpuarch=NA, minosversion=NA, queueforserver=FALSE, hostnode=NA, forcehost=FALSE, ramrequired=10, jobname=NA, jagspath='/usr/local/bin/jags', mgridpath=system.file("xgrid", "mgrid.sh", package="runjags"), hostname=Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), password=Sys.getenv("XGRID_CONTROLLER_PASSWORD"), ...){
	
	if(max.threads < 1) stop("The maximum number of threads must be greater than 0")
		
	passed <- c(list(...), model=model)
	if(any(names(passed)=="method") || any(names(passed)=="method.options" )) stop("Cannot specify 'method' or 'method.options' arguments to xgrid functions")
	if(any(names(passed)=="tempdir")) stop("Cannot specify 'tempdir' argument to xgrid.submit functions")
	if(any(names(passed)=="keep.jags.files")) stop("Cannot specify 'keep.jags.files' argument to xgrid.submit functions (this argument is supplied to xgrid.results.jags function instead)")
	
	# Translate monitor.pd/popt/pd.i/deviance (legacy code):
	monitor.deviance <- any(names(passed)=="monitor.deviance") && passed$monitor.deviance
	monitor.pd <- any(names(passed)=="monitor.pd") && passed$monitor.pd
	monitor.popt <- any(names(passed)=="monitor.popt") && passed$monitor.popt
	monitor.pd.i <- any(names(passed)=="monitor.pd.i") && passed$monitor.pd.i
	if(any(names(passed)=="check.conv")){
		warning("Use of the 'check.conv' argument is deprecated - use 'summarise' to achieve some of the same function", call.=FALSE)
		if(!any(names(passed)=="summarise")) passed <- c(passed, list(summarise=passed$check.conv))
	}
	
	if(any(names(passed)=="monitor")){
		passed$monitor <- c(passed$monitor, if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")
	}else{
		if(any(c(monitor.deviance,monitor.pd,monitor.popt,monitor.pd.i))){
			passed <- c(passed, list(c(if(monitor.deviance) "deviance", if(monitor.pd) "pd", if(monitor.popt) "popt", if(monitor.pd.i) "pd.i")))
			warning("Use of the 'monitor.deviance', 'monitor.pd', 'monitor.popt' and 'monitor.pdi' arguments is deprecated - add the special variables 'dic', 'deviance', 'pd', 'popt' or 'pd.i' to the 'monitor' argument instead", call.=FALSE)
		}
	}
	
	# Get default argument values from specified functions for anything not given to this function:
 	arglist <- getargs(c('run.jags','add.summary'), passed)	
	
 	args <- arglist[which(names(arglist) %in% names(formals(setup.jagsfile)))]
	args$method <- if(any(.packages()=="rjags")) "rjags" else "simple"
	if(!any(names(args)=='data'))
		args$data <- parent.frame()
	if(!any(names(args)=='inits'))
		args$inits <- parent.frame()
 	obj <- do.call("setup.jagsfile", args=args)
	
	arglist <- arglist[names(arglist)!="method" & names(arglist)!="method.options" & names(arglist)!="tempdir" & names(arglist)!="keep.jags.files"]
	
	# Get xgrid method options:
	method.options <- setup.xgrid(separate=max.threads>1,JAGSversion=JAGSversion, Rversion="", packages=list(), artfun=function() writeLines("1"), email=email, profiling=profiling, cpuarch=cpuarch, minosversion=minosversion, queueforserver=queueforserver, hostnode=hostnode, forcehost=forcehost, ramrequired=ramrequired, jobname=jobname, cleanup=FALSE, showprofiles=FALSE, jagspath=jagspath, Rpath='/usr/bin/R', Rbuild='64', max.filesize="1GB", mgridpath=mgridpath, hostname=hostname, password=password, submitandstop=TRUE, jagsrun=TRUE)
	method.options$max.threads <- max.threads
	method.options <- method.options[c("jobname", "command", "customart", "jagspath", "submitandstop", "max.threads", "mgridpath", "hostname", "password")]
	
	if(max.threads>1 && any(obj$monitor %in% c("dic", "pd", "popt", "pd.i"))){
		warning("Setting max.threads to 1 to allow calculation of dic/pd/popt/pd.i")
		max.threads <- 1
	}
	method.options$max.threads <- max.threads
	
	
	# Slight hack - we're not calling actually run.jags so add combine=FALSE to the argument list:
	arglist <- c(arglist, list(runjags.object=obj, combine=FALSE, method="xgrid", method.options=method.options, keep.jags.files=TRUE, tempdir=FALSE))
	
 	args <- arglist[which(names(arglist) %in% c(names(formals(extend.jags)), names(formals(add.summary))))]
	if(!any(names(arglist)=="burnin")) arglist$burnin <- 5000
	
 	res <- do.call("extend.jags", args=args)
	
	warning("As of runjags version 2.0, the xgrid.xxx.jags functions are untested. If you are still using them, please email the package maintainer - otherwise they are likely to be removed from a future release.")
	return(res)
}


xgrid.results.jags <- function(background.runjags.object, wait=TRUE, cleanup=TRUE){
	
	if(class(background.runjags.object)!="runjagsbginfo") stop("An object produced by a background runjags method must be supplied (see the manual page for more details)")
	if(background.runjags.object$method!="xgrid"){
		stop("This JAGS process was not started using an xgrid method")
	}
	
	jobinfo <- list(directory=background.runjags.object$directory, jobname=background.runjags.object$jobname, jobid=background.runjags.object$jobid)
	
	swcat("Retrieving xgrid JAGS results...\n\n")
	
	# First go to xgrid.retrieve which checks the jobname exists and gets files back from xgrid (silent is always FALSE as if silent.jags was specified, then jags output is binned):
	output <- xgrid.retrieve(jobinfo=jobinfo, wait=wait, silent=FALSE, cleanup=cleanup, partialretrieve=FALSE, jags=TRUE)

	if(!output$done & wait) stop("There was an error running the model")
	if(!output$done & !wait) stop("The model has not yet finished running")
	
	# tempdir is always FALSE, keep.jags.files can't be specified to xgrid.submit functions:
	tempdir <- FALSE
	keep.jags.files <- !cleanup

	warning("As of runjags version 2.0, the xgrid.xxx.jags functions are untested. If you are still using them, please email the package maintainer - otherwise they are likely to be removed from a future release.")




	# Then copy/paste from extend.jags:
	
	echo <- FALSE
	starttime <- Sys.time()
	
	# I switched the argument names:
	sub.chains <- NA    # either FALSE (stop if any crashed), TRUE (try to recover), NA (use runjags.getOption), or a numeric vector of which chains to read in
	sub.samples <- NA	# either NA (use sample size) or a number of iterations to thin to
	if(is.na(sub.samples)) sub.samples <- FALSE
	
	if(class(background.runjags.object)=='character' && file.exists(file.path(background.runjags.object, 'jagsinfo.Rsave'))){
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
	
	# This was done in extend.jags before, but now we may override the original request for combine so do it here instead:
	if(background.runjags.object$combine) background.runjags.object$burnin <- 0 else background.runjags.object$burnin <- background.runjags.object$burnin + background.runjags.object$adapt
	
	if(identical(sub.samples, TRUE) || identical(sub.samples, NA)) stop('Value for "sub.samples" must be either FALSE or the number of iterations to return')
	
	if(identical(sub.chains, NA)){
		if(!background.runjags.object$combine && runjags.getOption("partial.import")) sub.chains <- TRUE else sub.chains <- FALSE
	}

	if(!identical(sub.chains, FALSE) && background.runjags.object$combine) stop('Cannot return partially completed jobs when combining MCMC objects - try again with combine=FALSE')
	
	
	# Make name easier to type:
	runjags.object <- background.runjags.object
	
	summaryargs <- getsummaryargs(runjags.object$summary.pars)		
	
	if(!is.na(keep.jags.files)) runjags.object$keep.jags.files <- keep.jags.files
	if(!identical(read.monitor, NA)){
		read.monitor <- checkvalidmonitorname(read.monitor)
		# Switch things in monitor not in read.monitor to noread.monitor:
		tomove <- runjags.object$monitor[! runjags.object$monitor %in% read.monitor]
		# Preserve these (not really monitors anyway):
		tokeep <- runjags.object$monitor[runjags.object$monitor %in% c('pd','pd.i','popt','dic')]
		runjags.object$monitor <- c(read.monitor, tokeep)
		runjags.object$noread.monitor <- c(runjags.object$noread.monitor, tomove)
		# If read.monitor is a subset of an array, the whole array will be moved to noread and the subset kept in monitor (but with indices expanded)
	}
	
	# Call runjags.readin and then deal with copying files etc:
	
	if(length(runjags.object$noread.monitor)>0){
		read.monitor <- runjags.object$monitor[!runjags.object$monitor%in%c('pd','pd.i','popt','dic')]
	}else{
		read.monitor <- NA
	} 
	
	allok <- FALSE
	# This will be called if runjags.readin fails (allok=FALSE) or at the end of the function (allok=TRUE):
	on.exit({
		
		new.directory <- runjags.object$directory
		
		if(runjags.object$keep.jags.files){
			if(!allok) swcat('One or more simulations failed - you may be able to retrieve any successful simulations using:\nresults.jags("', new.directory, '", recover.chains=TRUE)\nSee the help file for that function for possible options.\n', sep='')
			if(!new.directory %in% runjagsprivate$simfolders) runjagsprivate$simfolders <- c(runjagsprivate$simfolders, new.directory)
		}else{
			if(!allok && runjags.getOption('keep.crashed.files') && new.directory!="Directory not writable"){
				swcat('One or more simulations failed so have not been deleted - you may be able to retrieve any successful simulations using:\nresults.jags("', new.directory, '", recover.chains=TRUE)\nSee the help file for that function for possible options.\n', sep='')
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
	
	simtimetaken <- difftime(newoutput$end.time, runjags.object$startedon)
	
	
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
	if(runjags.object$combine){
		burnin <- runjags.object$oldburnin

		combinedoutput <- list(mcmc=combine.mcmc(list(runjags.object$mcmc, newoutput$mcmc), collapse.chains=FALSE), deviance.table=weightedaverage(runjags.object$deviance.table, newoutput$deviance.table, niter(runjags.object$mcmc), niter(newoutput$mcmc)), deviance.sum=weightedaverage(runjags.object$deviance.sum, newoutput$deviance.sum, niter(runjags.object$mcmc), niter(newoutput$mcmc)), pd=if('pd'%in%runjags.object$monitor) combine.mcmc(list(runjags.object$pd, newoutput$pd), collapse.chains=FALSE) else NA, end.state=newoutput$end.state, samplers=newoutput$samplers)

	}else{
		combinedoutput <- list(mcmc=newoutput$mcmc, deviance.table=newoutput$deviance.sum, deviance.table=newoutput$deviance.sum, pd=newoutput$pd, end.state=newoutput$end.state, samplers=newoutput$samplers)
		# Never change thin - it messes everything up. Iteration labels will be slightly wrong, that's all.
		# runjags.object$thin <- thin
		# This means if extended, the thin will be what we had here
	}

	# Save some RAM:
	rm(newoutput)
	gcinfo <- gc(FALSE)

	# Takes into account (1) the time taken to import, (2) the time taken for the previous sims if any, (3) the simulation run time NOT including any gap between finishing and importing
	timetaken <- (difftime(Sys.time(), starttime) + runjags.object$timetaken + simtimetaken)
	
	if(sample != niter(combinedoutput$mcmc)) warning(paste('Unexpected discrepancy in sample size: expected ', sample, ' iterations but there are actually ', niter(combinedoutput$mcmc), sep=''))
	
	# Call function to calculate summary statistics and plots etc:	
	combinedoutput <- makerunjagsobject(combinedoutput, summarise=runjags.object$summarise, summaryargs=summaryargs, burnin=burnin, sample=sample, thin=thin, model=runjags.object$model, data=runjags.object$data, monitor=runjags.object$monitor, noread.monitor=runjags.object$noread.monitor, modules=runjags.object$modules, factories=runjags.object$factories, response=runjags.object$response, residual=runjags.object$residual, fitted=runjags.object$fitted, method=runjags.object$method, method.options=runjags.object$method.options, timetaken=timetaken)

	stopifnot(class(combinedoutput$end.state)=='runjagsinits')
	
	swcat("Finished running the simulation\n")
	
	return(combinedoutput)

	
}


xgrid.run.JAGS <- xgrid.run.jags
xgrid.autorun.JAGS <- xgrid.autorun.jags
xgrid.extend.JAGS <- xgrid.extend.jags
xgrid.autoextend.JAGS <- xgrid.autoextend.jags
xgrid.submit.JAGS <- xgrid.submit.jags

xgrid.results.JAGS <- xgrid.results.jags


# Deprecated (but no warnings as they are just copied):
xgrid.run.jagsfile <- xgrid.run.jags
xgrid.autorun.jagsfile <- xgrid.autorun.jags
xgrid.submit.jagsfile <- xgrid.submit.jags
xgrid.run.JAGSfile <- xgrid.run.jagsfile
xgrid.autorun.JAGSfile <- xgrid.autorun.jagsfile
xgrid.submit.JAGSfile <- xgrid.submit.jagsfile
