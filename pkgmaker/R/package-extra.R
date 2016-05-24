# Package extra action registry
# 
# Author: renaud
###############################################################################

#' @include devutils.R 
NULL

.getExtraEnv <- function(package){
	if( missing(package) || is.null(package) ) where <- topns(FALSE)
	else if( isString(package) ) {
		package <- sub("^package:", "", package)
		if( package == 'R_GlobalEnv') where <- .GlobalEnv
		else where <- asNamespace(package)
	}
	else stop("Invalid argument `package`: must be missing or a package name.")
	where
}

# extra handler registry
extra_handlers <- setPackageRegistry('extra_handler', 'function' 
									, description = 'Handler functions for package-specific extra tasks'
									, entrydesc = 'extra handler')
							
# extra action registry
extra_actions <- registry()
extra_actions$set_field("key", type="character", is_key = TRUE, index_FUN = match_exact)
extra_actions$set_field("package", type="character", is_key = TRUE, index_FUN = match_exact)
extra_actions$set_field("handler", type='character', is_mandatory=TRUE, is_key=TRUE)
extra_actions$set_field("args", type='list', default=list())
extra_actions <- setPackageRegistry('extra_action', extra_actions
									, description = 'Handler functions for package-specific extra actions'
									, entrydesc = 'extra action')

#' Install/Run Extra Things After Standard Package Installation
#' 
#' @description
#' These functions define a framework to register actions for which default sets of arguments
#' can be defined when (lazy-)loading a package, and run later on, e.g., after the package 
#' is installed using dedicated commands.
#' 
#' \code{setPackageExtraHandler} defines main action handler functions, for which 
#' actions are defined as a set of arguments and registered using \code{setPackageExtra}. 
#'  
#' @param handler name of a handler, e.g, \code{'install'}.
#' It must be unique across all handlers registered by any other packages.  
#' @param fun handler function that will be called with the arguments registered
#' with \code{packageExtra(name, ...)}
#' @param package package name where to store/look for the internal registries.
#' End users should not need to use this argument.
#' 
#' @return the runner function associated with the newly registered handler,
#' as built by \code{packageExtraRunner}.  
#'  
#' @rdname packageExtra
#' @export
setPackageExtraHandler <- function(handler, fun, ...){
	
	# add entry to the registry
	setPackageRegistryEntry('extra_handler', handler, fun, ...)
	# build associated runner
	runner <- packageExtraRunner(handler)
}

#' \code{packageExtraHandler} retrieves a given handler from the registry. 
#' 
#' @param ... extra arguments passed to internal function calls.
#' In \code{packageExtraHandler}, these are passed to \code{\link{pkgreg_fetch}}.
#' 
#' In \code{setPackageExtra}, these define default arguments for the handler function. 
#' These are overwritten by arguments in the call to runner function if any.
#'  
#' @rdname packageExtra
#' @export 
packageExtraHandler <- function(handler=NULL, ...){
	# load handler from registry
	pkgreg_fetch('extra_handler', key=handler, ...)
}
#' \code{setPackageExtra} registers extra actions for a given handler.
#' 
#' For example, calling \code{setPackageExtra('install', pkgs='non_CRAN_pkg', repos='http://non-standard-repo')}
#' in a source file of package 'myPkg' registers the call 
#' \code{install.packages('non_CRAN_pkg', repos='http://non-standard-repo', ...)}
#' in a registry internal to the package. 
#' All calls to \code{setPackageExtra('install', ...)} can then be run by the user, as
#' a post installation step via \code{install.extrapackages('myPkg', ..)}.
#' 
#' @param extra name of the extra action.
#' @param .wrap logical that indicates if a function that runs the extra action should
#' be returned or only the default arguments
#' 
#' @rdname packageExtra
#' @export
setPackageExtra <- function(handler, extra, ...){
	
	# check that a handler is defined in the registry
	fhandler <- packageExtraHandler(handler, exact=TRUE, error=FALSE)
	if( is.null(fhandler) ){
		handlers <- packageExtraHandler()
		stop("Could not register action '", extra, "': handler '", handler, "' is not defined"
				, if( length(handlers) ){
					str_c(".\n  Available handlers are: ", str_out(handlers, Inf))
				} else " [handler registry is empty]." )
	}
	args <- list(...)
	pkg <- packageName(topenv(parent.frame()), .Global=TRUE)
	setPackageRegistryEntry('extra_action', key=extra, handler=handler, args=args
							, package = pkg 
							, msg=str_c(" for handler '", handler, "'"))
}


.wrapExtra <- function(fhandler, args=list()){
	
	# define wrapper function
	f <- function(...){
		cl <- match.call()
		cl[[1L]] <- as.name('fhandler')
		# add default arguments
		lapply(names(args), function(a){
			if( !a %in% names(cl) )
				cl[[a]] <<- as.name(substitute(a, list(a=a)))
		})
		eval(cl)
	}
	# set registered arguments as default arguments
	formals(f) <- c(args, formals(f))
	f
}
#' \code{packageExtra} retrieve a given extra action, either as its registry entry,
#' or as a function that would perform the given action.
#' @rdname packageExtra
#' @export
packageExtra <- function(handler=NULL, extra=NULL, package=NULL, .wrap=FALSE){
	
	# load extra registry
	extras <- pkgreg_fetch('extra_action', key=extra, handler=handler, package=package
						, exact=TRUE, all=!.wrap)
	
	# return whole registry if no other argument is provided
	if( missing(handler) || is.null(extra) || !.wrap ) return( extras )
		
	args <- extras$args
	fhandler <- packageExtraHandler(handler, package='pkgmaker')
	if( is.null(fhandler) ){
		handlers <- packageExtraHandler(package='pkgmaker')
		stop("Could not find action handler '", handler, "' in pkgmaker global handler registry.\n"
				, "  Available handlers are: ", str_out(handlers, Inf))
	}
	# define wrapper function
	.wrapExtra(fhandler, args)		
}
#' \code{packageExtraRunner} defines a function to run all or some of the actions registered 
#' for a given handler in a given package.
#' For example, the function \code{install.extrapackages} is the runner defined for the extra handler \code{'install'} 
#' via \code{packageExtraRunner('install')}.
#' 
#' @param .verbose logical that indicates if verbose messages about the extra actions being
#' run should be displayed.
#' 
#' @rdname packageExtra
#' @export
packageExtraRunner <- function(handler){

	.handler <- handler
	function(package, extra=NULL, handler=NULL, ..., .verbose=getOption('verbose')){
		
		if( missing(handler) ) handler <- .handler
		.local <- function(p, ...){
			# load list of extras
			extras <- packageExtra(handler=handler, extra=extra, package=p)
			# execute extras
			sapply(extras, 
				function(def, ...){
					e <- def$key
					h <- def$handler
					f <- packageExtra(handler=h, extra=e, package=p, .wrap=TRUE)
					if( .verbose ){
						message("# Running extra action '", h, ':', e, "' ...")
						message("# Action: ", str_fun(f))
						on.exit( message("# ERROR [", e, "]\n") )
					}
					res <- f(...)
					if( .verbose ){
						on.exit()
						message("# OK [", e, "]\n")
					}
					res
				}
			, ...)
		}
		invisible(sapply(package, .local, ...))
	}
}

#' \code{install.extrapackages} runs all extra actions registered for a given package.
#' 
#' @rdname packageExtra
#' @export
install.extras <- packageExtraRunner(NULL)
#' \code{install.extrapackages} install sets of packages that can enhance a 
#' package, but may not be available from CRAN.
#' 
#' \code{install.extrapackages} is defined as the extra handler for 
#' the extra action handler \code{'install.packages'}.
#' All arguments in \code{...} are passed to \code{\link{install.packages}}.
#' By default, packages that are already installed are not re-installed.
#' An extra argument \code{force} allows to force their installation.
#' The packages are loaded if their installation is successful. 
#' 
#' @rdname packageExtra
#' @export
install.extrapackages <- setPackageExtraHandler('install.packages', 
	function(pkgs, ..., force=FALSE){
		res <- sapply(pkgs, function(pkg, ...){
			if( force || !require.quiet(pkg, character.only=TRUE) ){
				install.packages(pkg, ...)
				require(pkg, character.only=TRUE)
			}else message("Loaded extra package: ", pkg)
		}, ...)
	}
)
