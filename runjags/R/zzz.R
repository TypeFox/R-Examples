.onLoad <- function(libname, pkgname){
	
	runjagsprivate$runjagsversion <- utils::packageDescription(pkgname, fields='Version')
	
	# Get and save the library location, getting rid of any trailing / caused by r_arch being empty:
	modloc <- gsub('/$','', file.path(libname, pkgname, 'libs', if(.Platform$r_arch!="") .Platform$r_arch else ""))
	if(!file.exists(file.path(modloc, paste('runjags', .Platform$dynlib.ext, sep=''))))
		modloc <- ''
	runjagsprivate$modulelocation <- modloc
	
	setopts <- mget('.runjags.options', envir=.GlobalEnv, ifnotfound=list(.runjags.options=NULL))[[1]]
	if(!is.null(setopts)){
		if(!is.list(setopts)){
			warning('Ignoring invalid (non-list) specification for .runjags.options on loading the runjags package', call.=FALSE)
		}else{
			newopts <- do.call('runjags.options', args=setopts)
		}
	}
	
	# Catch GUI and Rstudio specific settings: - not for Rgui (windows) though
	if(.Platform$GUI%in%c("AQUA") || Sys.getenv('RSTUDIO')!="")
		runjags.options(new.windows=FALSE)
	
	# To ensure that cleanup.jags is run when R is quit:
	reg.finalizer(runjagsprivate, .onDetach, onexit = TRUE)
}

.onAttach <- function(libname, pkgname){
	
	# This will be run after load if the package is attached:
	setopts <- mget('.runjags.options', envir=.GlobalEnv, ifnotfound=list(.runjags.options=NULL))[[1]]
	if(!is.null(setopts) && !runjags.getOption('silent.runjags')){
		packageStartupMessage(paste('Attaching runjags (version ', runjagsprivate$runjagsversion, ') and setting user-specified options', sep=''))
	}
}

.onDetach <- function(libpath){
	# Just in case it is not always safe to try and access an element of an env that is in the process of being deleted (when R quits):
	if(!is.null(runjagsprivate$dynlibname)) dynunloadmodule()
	all.folders <- try(runjags.getOption('full.cleanup'), silent=TRUE)
	if(class(all.folders)=='try-error') all.folders <- FALSE
	try(cleanup.jags(all.folders=all.folders, silent=TRUE), silent=TRUE)
}
