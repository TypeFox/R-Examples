load.runjagsmodule <- function(fail=TRUE, silent=FALSE){
	
	intfun <- function(){
		# Make sure this is not the module-less version from sourceforge:
		if(runjagsprivate$modulelocation=='')
			return('The internal runjags module is not installed - please reinstall the full version of the package from CRAN, or alternatively you can download a standalone version of the JAGS module from the sourceforge page at http://sourceforge.net/projects/runjags/')
	
		# Also makes sure JAGS is installed:
		if(!loadandcheckrjags(FALSE)) 
			return("The rjags package is required to use the internal runjags module - alternatively you can download a standalone version of the JAGS module from the sourceforge page at http://sourceforge.net/projects/runjags/")
	
		# Check the JAGS major version is as expected:	
		if(packageVersion('rjags')$major < runjagsprivate$minjagsmajor)
			return(paste('JAGS version ', runjagsprivate$minjagsmajor, '.x.x to ', runjagsprivate$maxjagsmajor, '.x.x is required for this version of the runjags module - please update JAGS and rjags',sep=''))
		if(packageVersion('rjags')$major > runjagsprivate$maxjagsmajor)
			return(paste('This version of the runjags module was designed for JAGS version ', runjagsprivate$minjagsmajor, '.x.x to ', runjagsprivate$maxjagsmajor, '.x.x - please update the runjags package', sep=''))
	
		success <- try(rjags::load.module('runjags',runjagsprivate$modulelocation))
	
		if(inherits(success, 'try-error')){
			
			rvers <- paste('version ', R.version$major, sep='')
			if(grepl('mac.binary', .Platform$pkgType, fixed=TRUE)){
				# A specific error may be because of SL vs Mavericks version on OS X for JAGS version 3.4:
				mavericks <- grepl('mavericks', .Platform$pkgType)
				if(mavericks)
					rvers <- paste(rvers, ' - Mavericks', sep='')
				else
					rvers <- paste(rvers, ' - Snow Leopard', sep='')
			}
		
			return(paste("The internal runjags module could not be loaded - perhaps the package was not built using the same versions of R [", rvers, "] and JAGS [version ", testjags(silent=TRUE)$JAGS.version, "] as available on this system?", sep=''))

		}
		return(TRUE)
	}
	
	retval <- intfun()
	
	if(retval==TRUE){
		invisible(TRUE)
	}else{	
		if(fail)
			stop(retval)
	
		if(!silent)
			swcat(retval,'\n',sep='')
	
		invisible(FALSE)
	}
}

unload.runjagsmodule <- function(){
	
	if(!loadandcheckrjags(FALSE))
		stop("The rjags package is required to use the internal runjags module - alternatively you can download a standalone version of the JAGS module from the sourceforge page at http://sourceforge.net/projects/runjags/")
		
	suppressWarnings(success <- try(rjags::unload.module('runjags')))
	
	if(inherits(success, 'try-error')){
		warning("There was a problem unloading the internal runjags module - if you installed this package from CRAN, please file a bug report to the package author")
		invisible(FALSE)
	}else{
		invisible(TRUE)
	}
}

load.runJAGSmodule <- load.runjagsmodule
unload.runJAGSmodule <- unload.runjagsmodule


# These utility functions are NOT exported, and are primarily used for unit testing.
# Availability and/or operation of these functions may change without warning.

dynloadmodule <- function(){
	
	# Sets environmental variables we need for Windows:
	if(.Platform$OS.type=='windows'){
		if(!loadandcheckrjags(FALSE))
			stop('The rjags package is required to load the internal dynlib')
	}
	
	if(runjagsprivate$modulelocation==''){
		warning('The runjags module has not been installed with this version of the package - try again using the CRAN binary')
		invisible(FALSE)
	}
	
	# Check the JAGS major version is as expected:
	if(testjags(silent=TRUE)$JAGS.major < runjagsprivate$minjagsmajor)
		return(paste('JAGS version ', runjagsprivate$minjagsmajor, '.x.x to ', runjagsprivate$maxjagsmajor, '.x.x is required for this version of the runjags module - please update JAGS and rjags',sep=''))
	if(testjags(silent=TRUE)$JAGS.major > runjagsprivate$maxjagsmajor)
		return(paste('This version of the runjags module was designed for JAGS version ', runjagsprivate$minjagsmajor, '.x.x to ', runjagsprivate$maxjagsmajor, '.x.x - please update the runjags package', sep=''))

	# Find and load the runjags shared library (only required for these tests and using the rjags call 'load.modue()' so NOT loaded at runtime):
	slibpath <- file.path(runjagsprivate$modulelocation, paste('runjags', .Platform$dynlib.ext, sep=''))
	swcat("Loading shared library from:  ", slibpath, "\n", sep="")
	success <- try(dyn.load(slibpath))
	
	if(inherits(success, 'try-error')){
		
		rvers <- paste('version ', R.version$major, sep='')
		if(grepl('mac.binary', .Platform$pkgType, fixed=TRUE)){
			# A specific error may be because of SL vs Mavericks version on OS X for JAGS version 3.4:
			mavericks <- grepl('mavericks', .Platform$pkgType)
			if(mavericks)
				rvers <- paste(rvers, ' - Mavericks', sep='')
			else
				rvers <- paste(rvers, ' - Snow Leopard', sep='')
		}
	
		return(paste("The runjags dynlib could not be loaded - perhaps the package was not built using the same versions of R [", rvers, "] and JAGS [version ", testjags(silent=TRUE)$JAGS.version, "] as available on this system?", sep=''))

	}
	
	runjagsprivate$dynlibname <- success
	invisible(TRUE)

}

dynunloadmodule <- function(){
	
	if(is.null(runjagsprivate$dynlibname)){
		warning('Unable to load the dynlib as it has not been loaded')
		invisible(FALSE)
	}
	# Find and unload the runjags shared library (only required for these tests and using the rjags call 'load.modue()' so NOT loaded at runtime):
	slibpath <- system.file("libs", paste(.Platform$r_arch, if(.Platform$r_arch!="") "/" else "", if(.Platform$OS.type=="unix") "runjags.so" else "runjags.dll", sep=""), package="runjags")
	swcat("Unloading shared library from:  ", slibpath, "\n", sep="")
	success <- try(dyn.unload(slibpath))
	if(inherits(success, 'try-error'))
		stop("The internal dynlib could not be unloaded - if you installed this package from CRAN, please file a bug report to the package author")

	runjagsprivate$dynlibname <- NULL
	invisible(TRUE)
}

userunjagsmodule <- function(distribution, funtype, parameters, x, uselog=FALSE, lower=TRUE){
	
	if(is.null(runjagsprivate$dynlibname)){
		stopifnot(dynloadmodule())
	}
	
	if(class(distribution)!="character") stop("The distribution type must be one of par1, par2, par3, par4, lomax, mouch, genpar or halfcauchy")
	disttype <- switch(distribution, par1=1, par2=2, par3=3, par4=4, lomax=5, mouch=6, genpar=7, hcauchy=8, halfcauchy=8, 0)
	if(disttype==0) stop("The distribution type must be one of par1, par2, par3, par4, lomax, mouch, genpar or halfcauchy")

	if(class(funtype)!="character") stop("The function type must be one of d, p and q")
	dpqr <- switch(funtype, d=1, p=2, q=3, r=4, 0)
	if(dpqr==0) stop("The function type must be one of d, p and q")
	if(dpqr==4) stop("The function type must be one of d, p and q - r is not available")

	npars <- length(parameters)
	N <- length(x)
	values <- numeric(N)

	# Having problems with bools on sparc, so use ints here:
	output <- .C('testrunjags',PACKAGE='runjags',disttype=as.integer(disttype),dpqr=as.integer(dpqr),uselog=as.integer(uselog),lower=as.integer(lower),N=as.integer(N), x=as.double(x),npars=as.integer(npars),parameters=as.double(parameters),values=as.double(values),status=integer(1))

	if(output$status==1) stop("Incorrect number of parameters provided")
	if(output$status==2) stop("Unrecognised distribution type")
	if(output$status==3) stop("Invalid parameter values provided")
	if(output$status==4) stop("Function type not valid")
	if(output$status==5) stop("Function type not valid")

	return(output$values)

}
