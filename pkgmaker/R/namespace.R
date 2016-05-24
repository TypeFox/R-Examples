# Namespace related functions
# 
# Author: Renaud Gaujoux
# Creation: 30 Apr 2012
###############################################################################


is_funcall <- function(fun){

	n <- sys.nframe()
	i <- 1
	dg <- digest(fun)
	while( i <= n ){
		f <- sys.function(i)
		ca <- sys.call(i)
#		cat(digest(f), dg, getPackageName(environment(f), FALSE), "\n")
		if( digest(f) == dg ) return(i)
		i <- i + 1
	}
	FALSE
}

is_pkgcall <- function(pkg){
	
	pkg %in% pkg_calls()
	
}

pkg_calls <- function(){
	n <- sys.nframe() - 1
	i <- 1
	res <- character()
	while( i <= n ){
		f <- sys.function(i)
		e <- environment(f)
		if( !is.null(e) ){
			pkg <- methods::getPackageName(e, create=FALSE)
			if( pkg != '' ) res <- c(res, pkg)
		}
		i <- i + 1
	}
	res
}

#' Namespace Development Functions
#' 
#' \code{getLoadingNamespace} returns information about the loading namespace.
#' It is a wrapper to \code{\link{loadingNamespaceInfo}}, that does not throw 
#' an error.
#' 
#' @param env logical that indicates that the namespace's environment (i.e. the 
#' namespace itself) should be returned.
#' @param info logical that indicates that the complete information list should 
#' be returned
#' 
#' @return the name of the loading namespace if \code{env} and \code{info} are 
#' \code{FALSE}, an environment if \code{env=TRUE}, a list with elements 
#' \code{pkgname} and \code{libname} if \code{info=TRUE}. 
#' 
#' @rdname namespace
#' @export
#' 
getLoadingNamespace <- function(env=FALSE, info=FALSE, nodev=FALSE){
	is.loading <- try(nsInfo <- loadingNamespaceInfo(), silent=TRUE)
	if( !is(is.loading, 'try-error') ){
		if( env ) asNamespace(as.name(nsInfo$pkgname))
		else if( info ){
			nsInfo$path <- file.path(nsInfo$libname, nsInfo$pkgname)
			nsInfo 
		}else nsInfo$pkgname
		
	}else if( !nodev ){ # devtools namespaces are allowed
		if( is_pkgcall('devtools') && (i <- is_funcall(devtools::load_all)) ){
			# find out the package that is currently being loaded by load_all
			e <- sys.frame(i)
			pkg <- e$pkg
			# extract namespace
			if( env ) asNamespace(pkg$package)
			else if( info ){
				list(
						pkgname = pkg$package
						, path = pkg$path
						, libname = dirname(pkg$path)
				)
			}else pkg$package
		}
	}
	else NULL
}

#' Tests if a namespace is being loaded.
#' 
#' @param ns the name of a namespace or a namespace whose loading state is tested.
#' If missing \code{isLoadingNamespace} test if any namespace is being loaded.
#' @param nodev logical that indicates if loading devtools namespace should 
#' be discarded.
#' 
#' @rdname namespace
#' @export
isLoadingNamespace <- function(ns, nodev=FALSE){
	
	if( missing(ns) ) !is.null(getLoadingNamespace(nodev=nodev))
	else{
		nspkg <- getLoadingNamespace(nodev=nodev, env=is.environment(ns))
		if( is.null(nspkg) ) FALSE
		else identical(nspkg, ns)
	}
}

#' \code{isNamespaceLoaded} tests if a given namespace is loaded, without loading it, 
#' contrary to \code{\link{isNamespace}}.
#' 
#' @rdname namespace
#' @export
isNamespaceLoaded <- function(ns){
	if( is.environment(ns) ){
		if( !isNamespace(ns) ) return(FALSE)
		else ns <- getPackageName(ns)
	}
	if( isString(ns) ) ns %in% loadedNamespaces()
	else stop("Invalid argument `ns`: only support strings and environments.")
}

#' \code{isDevNamespace} tests the -- current -- namespace is a devtools namespace.
#' 
#' @rdname namespace
#' @export
isDevNamespace <- function(ns){
	if( missing(ns) ){
		e <- parent.frame()
		ns <- methods::getPackageName(topenv(e))
	}
	
	# cannot be true if the namespace is not loaded
	if( !isNamespaceLoaded(ns) ) return( FALSE )
	# get the namespace environment
	if( isString(ns) ) ns <- asNamespace(ns)
	# check for the presence of a .__DEVTOOLS__ object
	exists('.__DEVTOOLS__', where=ns)
	
}

#' Dynamically adds exported objects into the loading namespace.   
#' 
#' @param x character vector containing the names of R objects to export in the 
#' loading namespace.
#' 
#' @rdname namespace
#' @export
addNamespaceExport <- function(x){
	ns <- pkgmaker::getLoadingNamespace(env=TRUE)
	if( !is.null(ns) ){
		namespaceExport(ns, x)
	}
}

#' \code{ns_get} gets an object from a given namespace.
#' @rdname namespace
#' @export
ns_get <- function(x, ns){
    if( !isNamespace(ns) ) ns <- asNamespace(ns)
    get(x, ns)
}
