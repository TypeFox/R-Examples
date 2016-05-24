# Package hooks
# 
# Author: renaud
# Creation: 26 Jun 2012
###############################################################################

#' @include utils.R
#' @include devutils.R
#' @import stats
#' @import methods
NULL

#' Default Load/Unload Functions
#' 
#' @inheritParams base::.onLoad
#' @inheritParams base::library.dynam
#' 
#' @export
#' @rdname load
#' 
#' @examples
#' 
#' # in a package namespace:
#' .onLoad <- function(libname=NULL, pkgname){
#' 
#' 	pkgmaker::onLoad(libname, pkgname)
#' 
#' }
onLoad <- function(libname=NULL, pkgname, chname=packageName()){
	
	# load compiled library normally or in devmode
	if( !is.null(libname) ){
		if( file.exists(packagePath('libs')) ){
			sapply(chname, library.dynam, package=pkgname, lib.loc=libname)
		}
	}else{
		compile_src() # compile source files and load
	}
		
}

#' @inheritParams base::.onUnload
#' @export
#' @rdname load
#' 
#' @examples
#' 
#' # in a package namespace:
#' .onUnload <- function(libpath){
#' 
#' 	pkgmaker::onUnload(libpath)
#' 
#' }
onUnload <- function(libpath) {
	
	# unload compiled library normally or in devmode
	dlls <- base::getLoadedDLLs()
	pname <- packageName()
	if ( pname %in%  names(dlls) ){
		if( !missing(libpath) )	library.dynam.unload(pname, libpath)
		else dyn.unload(dlls[[pname]][['path']])
	}
	
}


#' Postponing Actions
#' 
#' This function implement a mechanism to postpone actions, which can be executed
#' at a later stage.
#' This is useful when developing packages, where actions that need to be run in the 
#' \code{link{.onLoad}} function but can be defined close to their context.
#' 
#' @param expr expression that define the action to postpone.
#' Currently only functions are supported.
#' @param key identifier for this specific action.
#' It should be unique across the postponed actions from the same group. 
#' @param group optional parent action group.
#' This enables to define meaningful sets of actions that can be run all at once.
#' @param envir environment in which the action should be executed.
#' Currently not used.
#' @param verbose logical that toggles verbose messages.
#'
#' @import digest
#' @export
#' 
#' @examples
#' opt <- options(verbose=2)
#' 
#' # define actions
#' postponeAction(function(){print(10)}, "print")
#' postponeAction(function(){print(1:10)}, "more")
#' postponeAction()
#' # execute actions
#' runPostponedAction()
#' runPostponedAction()
#' 
#' # restore options
#' options(opt)
#' 
postponeAction <- function(expr, key=digest(tempfile()), group=NULL, envir=topns(strict=FALSE), verbose=getOption('verbose')){
	
	# do not do anything if already running delayed actions
	if( isRunningPostponedAction() ) return()
	
	ns <- topns(strict=FALSE)
	taskObj <- simpleRegistry('.__delayedTasks__', envir=ns)
	if( !missing(expr) ){
		if( missing(key) ){
			stop("Missing required argument `key` for registering/cancelling delayed action.")
		}
		# add group prefix
		if( !is.null(group) )
			key <- str_c(group, '::', key)
		#qe <- if( !is.language(expr) ) substitute(expr) else expr
		qe <- expr
		if( verbose ){
			if( !is.null(qe) ) message("# Postponing action '", key, "'")
			else{
				message("# Cancelling postponed action '", key, "'")
			}
		}
		taskObj$set(key, list(action=qe, envir=envir))
	}else{
		taskObj$names()
	}
}

#' @rdname postponeAction
#' @export
runPostponedAction <- function(group=NULL, verbose=getOption('verbose')){
	
	ns <- topns(strict=FALSE)
	taskObj <- simpleRegistry('.__delayedTasks__', envir=ns)
	
	if( verbose ){
		message("# Executing postponed "
				, if( !is.null(group) ) paste("'", group, "' ", sep='')
				, "action(s) in package '"
				, packageName(ns, .Global=TRUE), "' ... "
				, appendLF = FALSE)
	}
	# set up running flag
	isRunningPostponedAction(TRUE)
	on.exit(isRunningPostponedAction(FALSE))
	#
	# execute actions
	t <- taskObj$names()
	if( !is.null(group) ) t <- grep(str_c("^", group), t, value=TRUE)
	if( verbose > 1 && length(t) ) message()
	sapply(t, function(x){
				act <- taskObj$get(x)
				if( verbose > 1 ){
					message("** Action '", x, "' [", packageName(act$envir, .Global=TRUE), ']')
				}
				act$action()
				taskObj$set(x, NULL)
				#eval(x$expr, x$envir)
			})
	if( verbose ) message('OK [', length(t), ']')
	invisible(length(t))
}

# Tells if one is executing deferred tasks via \code{onLoad}
isRunningPostponedAction <- sVariable(FALSE)

#' Simple Package Registry
#' 
#' @param name name of the registry object, with which it will
#' be assigned in \code{envir}.
#' @param envir environment where to store the registry object.
#' Defaults to the caller's top environment.
#' @param verbose logical that toggle a verbose message when 
#' the object is first created.
#' 
#' @export
simpleRegistry <- function(name, envir=topenv(parent.frame()), verbose=FALSE){
	
	# return stored instance if it exists
	if( exists(name, envir=envir) ){
		return( invisible(get(name, envir=envir)) )
	}
	
	if( verbose ) message("# Setup simple registry '", name, "' in ", packageName(envir, .Global=TRUE))
	.name <- name
	.envir <- envir
	.data <- list()
	
	.get <- function(x){
		if( .has(x) ){
			.data[[x]]
		}
	}
	.set <- function(x, value){
		if( is.null(value) ){
			if( .has(x) ){
				.data[[x]] <<- NULL
			}
		}else{
			.data[[x]] <<- value
		}
	}
	.has <- function(x){
		x %in% names(.data)
	}
	.cleanup <- function(){
		rm(.name, envir=.envir)
	}
	.names <- function(){
		names(.data)
	}
	.length <- function(){
		length(.data)
	}
	
	.obj <- list(get=.get, set=.set, has=.has
			, cleanup=.cleanup, names=.names
			, length = .length)
	
	# assign container object
	assign(.name, .obj, envir=.envir)
	#
	invisible(.obj)
}


#' Defunct Functions in pkgmaker
#' 
#' These functions have been defunct or superseded by other 
#' functions. 
#' 
#' @param ... extra arguments
#' 
#' @rdname pkgmaker-defunct
#' @name pkgmaker-defunct
NULL