# Registry utility functions
# 
# Author: Renaud Gaujoux
# Creation: 25 Apr 2012
###############################################################################

#' @include devutils.R
#' @include is.R
NULL

#' Fix Registry Access Functions 
#' 
#' Fixes some of the member functions for registry objects, as defined in the 
#' \pkg{registry} package.
#' The main fixed issue is due to the fact that key filtering does not return 
#' the correct entry, when an entry key is a prefix of another entry key,
#' even when passed the exact key. 
#' 
#' @param regobj registry object
#' 
#' @family registry-internals
#' @keywords internal
fix_registry <- function(regobj){
	
	# get private environment
	.REGENV <- environment(environment(regobj$n_of_entries)$f)
	
	# do not fix twice
	if( isFixed <- exists('.isFixed', .REGENV, inherits=FALSE) ){
		return(regobj)
	}
	
#	message("REGENV:\n", capture.output(print(ls(.REGENV, all=TRUE))))
#	message("env(delete_entry)\n", capture.output(print(ls(environment(environment(regobj$delete_entry)$f), all=TRUE))))
	# dummy variables for R CMD check
	PERMISSIONS <- 
	.get_entry_indices <-  
	.get_entry_names <- 
	SEALED_ENTRIES <-
	DATA <- 
	.delete_entry <- 
	NULL

	# .get_entries
	.get_entries <- get('.get_entries', .REGENV, inherits=FALSE)
	.get_fields <- get('.get_fields', .REGENV, inherits=FALSE)
	hook <- function(...){
		# remove fields that are not keys
		fld <- .get_fields()
		if( length(fld) ){
			keyfield <- names(fld[sapply(fld, function(x) isTRUE(x$is_key) )])
			index_fields <- list(...)
			if( !is.null(names(index_fields)) ){
				i <- match(keyfield, names(index_fields))
				index_fields <- index_fields[i[!is.na(i)]]
			}
		}
		do.call(.get_entries, index_fields)
	}
	assign('.get_entries', hook, .REGENV)
	#
	
	# fix bug in delete_entry
	hook <- function(...){
		key <- list(...)
		isString <- function(x) is.character(x) && length(x) == 1L
		if( length(key) == 1L && isString(key[[1L]]) ){
			
			errmsg <- paste0("Could not delete entry '", key[[1L]],"': ")
			if (!PERMISSIONS["delete_entries"]) 
				stop(errmsg, "deletion of entries denied due to restricted permissions.", call. = FALSE)
			entry_index <- .get_entry_indices(key)
			
			# fix: check for exact match (on full key)
			if( key[[1L]] %in% .get_entry_names() ){
				entry_index <- match(key[[1L]], .get_entry_names())
			}
			#end_fix
			if( !length(entry_index) ){
#				if( !quiet ) warning(errmsg, "not in registry.", immediate.=TRUE, call. = FALSE)
				return()
			}
			if (length(entry_index) != 1)
				stop(errmsg, "multiple matches.", call. = FALSE)
			if (entry_index %in% SEALED_ENTRIES) 
				stop(errmsg, "deletion of entry not allowed in sealed registry.", call. = FALSE)
			DATA[entry_index] <<- NULL
		} else .delete_entry(...)
	}
	environment(hook) <- .REGENV
	regobj$delete_entry <- hook
	#

	# fix bug in get_entry
	hook <- function(...){
		key <- list(...)
		isString <- function(x) is.character(x) && length(x) == 1L
		if( length(key) == 1L && isString(key[[1L]]) ){
			res <- .get_entries(...)
			if( key[[1L]] %in% names(res) )
				res[[key[[1L]]]]
			else res[[1L]] 
		}else .get_entries(...)[[1]]
	}
	environment(hook) <- .REGENV
	regobj$get_entry <- hook
	#

	# flag the registry as fixed
	assign('.isFixed', TRUE, .REGENV)
	# return fixed registry
	regobj
}


testreg <- function(){
	regobj <- registry()
	regobj$set_field("X", type = TRUE)
	regobj$set_field("Y", type = "character")
	regobj$set_field("index", type = "character", is_key = TRUE,
			index_FUN = match_partial_ignorecase)
	# fix
	regobj <- fix_registry(regobj)
	
	regobj$set_entry(X = TRUE, Y = "bla", index = "test")
	regobj$set_entry(X = TRUE, Y = "bloblo", index = "test2")
	regobj$set_entry(X = FALSE, Y = "foo", index = c("test", "bar"))
	
	regobj
}

#' Creates or Retrieves a Package Meta Registry
#' 
#' This function is used in \code{\link{setPackageRegistry}} and 
#' \code{\link{packageRegistry}} to create or query meta registries.
#' 
#' @keywords internal
.packageMetaRegistry <- function(package, quiet=FALSE, create=FALSE){
	
#	library(registry)
	metaregname <- '.packageRegistry'
	# get package environment
	e <- packageEnv(package)	
	# get namespace name
	nm <- packageName(e)
	
	# create registry environment if necessary
	if( !exists(metaregname, e, inherits=FALSE) ){
		if( !create ){
			if( quiet ) return(NULL)
			# throw error
			stop("Meta registry in package `", nm, "` does not exist.")
		}
#		if( !isLoadingNamespace(e) ){
#			stop("Can only create a package meta-registry when loading it"
#				," [loading namespace: ", if(!is.null(ns <- getLoadingNamespace()) ) ns else 'none', "].")
#		}
		message("Creating meta registry in package '", nm, "' ... ", appendLF=FALSE)
		# create registry object with special classes for the registry and entries
		meta <- registry(c(paste('package', nm, 'subregistry', sep='_'), 'package_subregistry')
				, c(paste('package', nm, 'metaregistry', sep='_'), 'package_metaregistry'))
		## set fields
		# access key
		meta$set_field("key", type="character", is_key = TRUE, index_FUN = match_exact)
		# sub-registry object
		meta$set_field("regobj", type="registry", is_mandatory = TRUE)
		# human readable description
		meta$set_field("description", type="character", is_mandatory = TRUE)
		# short object description
		meta$set_field("entrydesc", type="character", is_mandatory = TRUE)
		# parent package = owner of the primary registry
		meta$set_field("parent", type="character", default = '')
		# owner package (its value is forced)
		meta$set_field("package", type="character", default = nm, alternatives=nm)
		#
		
		# fix registry
		meta <- fix_registry(meta)
		# add package attribute
		attr(meta, 'package') <- nm
		# store within the calling package environment
		assign(metaregname,  meta, envir = e)
		message('OK')
	}
	# get package meta registry
	get(metaregname, envir=e, inherits = FALSE)
}

#' Package Registry
#' 
#' \code{packageRegistry} provides ways to create query package specific 
#' registries.
#' 
#' Package registries are organised in a meta-registry (a registry of registries) within a package's namespace. 
#' Each registry can be used to store sets of built-in or user-defined objects 
#' in an organised way, e.g. algorithms or datasets.
#' 
#' A package meta-registry is a \code{\link[registry:regobj]{registry}} object,
#' whose entries are \code{\link[registry:regobj]{registry}} objects themselves.
#' A sub-registry entry is defined by the following fields:
#' \describe{
#' \item{key}{The sub-registry's accession key/identifier (a character string).}
#' \item{regobj}{The sub-registry itself (a \code{registry} object)}
#' \item{description}{Human readable description of the purpose of the registry (a character string)}
#' \item{description}{Short human readable description of the type of entries (a character string)}
#' \item{package}{owner package, which is forced to be the package in which the meta registry
#' is defined.}
#' \item{parent}{The name of the package that holds the parent registry, which we 
#' call the primary package.
#' This field is non empty for cross-package registries, i.e. registries that 
#' derive from primary package's own registry.
#' Their entries are defined when (lazy-)loading the dependent package's namespace.}
#' }
#' 
#' Note that this function cannot be called from the global environment, but from 
#' a package namespace, e.g., when a package is lazy-loaded on installation or loaded
#' via the function \code{\link[devtools]{load_all}} from the \pkg{devtools} package.  
#'   
#' @param regname Name of a sub-registry, used as its identifier.
#' @param quiet a logical that indicates that one should return the (meta-)registry if it exists, 
#' or \code{NULL} otherwise, without throwing any error.
#' @param entry logical that indicates if the corresponding meta registry entry should 
#' be directly returned, without any other processing.
#' @param update logical that indicates if the package registry should be updated, by adding/removing 
#' entries from other loaded/unloaded packages. 
#' @param package package where to store or look for the registry.
#' @return a \code{\link[registry:regobj]{registry}} object or \code{NULL} (see argument 
#' \code{quiet}).
#' 
#' @import registry
#' @rdname registry
#' @export
packageRegistry <- function(regname=NULL, quiet=FALSE, entry=FALSE, update=!entry, package=topenv(parent.frame())){
	
#	library(registry)
	metaregname <- '.packageRegistry'
	name <- regname

	# get package environment
	e <- packageEnv(package)	
	# get namespace name
	nm <- packageName(e)
	
	# get package meta-registry
	pkgreg <- .packageMetaRegistry(package, quiet)
	
	# return meta registry if no name is specified
	if( is.null(name) )	return(pkgreg)
	else{
		if( is.null(pkgreg) ){
			if( quiet ) return(NULL)
			# throw error
			stop("Could not find registry '", name, "' in package `", nm, "`: meta registry does not exist.")	
		} 
		# retrieve sub-registry entry
		nm <- packageSlot(pkgreg) 
		reg <- regfetch(pkgreg, key=name, exact=TRUE, error=FALSE)
		
		if( is.null(reg) ){# not found
			if( quiet ) return(NULL)
			# throw error
			stop("Could not find registry `", name, "` in package `", nm, "`.")
		}else{
			# synchronise and substitute by primary sub-registry (if necessary)
			if( update ) reg <- .update_pkgreg(reg)
			# return plain registry entry if requested
			if( entry )	return(reg)
			# return sub-registry object
			reg$regobj
		}
	}
}

.update_pkgreg <- local({
	.cacheNS <- list()
	.cacheMD5 <- list()
	function(regentry){
		
		verbose <- getOption('verbose')
		
		# directly return entry if:
		# - one is loading the namespace of the package (primary or not)
		if( isLoadingNamespace(regentry$package) ) return(regentry) 		
		# - not a primary registry
		if( nchar(regentry$parent) > 0L ) return(regentry)
		
		primary <- regentry$package
		primaryreg <- regentry$regobj
		key <- regentry$key
		fullkey <- str_c(primary, '::', key)
		# sync if loaded packages changed
		hash <- digest(c(.cacheNS[[fullkey]], ns <- loadedNamespaces()))
#		print(ns)
#		print(.cacheNS)
		if( !identical(hash, .cacheMD5[[fullkey]]) ){
			if( verbose ) message('Updating registry ', fullkey, " ... ", appendLF=FALSE)
			# remove entries from unloaded packages
			if( length(.cacheNS[[fullkey]]) && length(notloaded <- setdiff(.cacheNS[[fullkey]], ns)) ){
				ndel <- sapply(notloaded, function(p){
					if( verbose > 1L ) message("\n Removing entries from package ", p, " ... ", appendLF=FALSE)
					e <- primaryreg$get_entry_names()
					n <- sapply(e, function(x){
						rec <- primaryreg$get_entry(x)
						if( rec$REGISTERINGpackage == p ){
							primaryreg$delete_entry(x)
							1L
						}else 0L
					})
					if( verbose > 1L ) message('OK [', sum(n), ']')
					sum(n)
				})
			}
			
			# list packages that have local versions of this registry
			reglist <- packageRegistries(fullkey)
#			print(reglist)
			pkgs <- names(reglist)
			# add entries from new packages into the primary registry
			if( length(miss <- setdiff(pkgs, .cacheNS[[fullkey]])) ){
				nadd <- sapply(miss, function(p){
					if( verbose > 1L ) message("\n Adding entries from package ", p, " ... ", appendLF=FALSE)
					reg <- packageRegistry(fullkey, package=p)
					e <- reg$get_entries()
					n <- sapply(e, function(x){
						# add entry if it does not exists already
						oldentry <- regfetch(primaryreg, KEYS=x, exact=TRUE, error=FALSE)
						if( is.null(oldentry) ){
							do.call(primaryreg$set_entry, x)
							1L
						}else 0L				
					})
					if( verbose > 1L ) message('OK [', sum(n), ']')
					sum(n)
				})
			}
			# store contributing packages and MD5 hash
			.cacheNS[[fullkey]] <<- pkgs
			.cacheMD5[[fullkey]] <<- digest(c(.cacheNS[[fullkey]], ns)) 
			if( verbose ) message('OK')
		}
		
		regentry
	}
})


#' \code{packageRegistries} lists registries from loaded packages.
#' 
#' @param primary logical that indicates if only primary registries 
#' should be listed.
#'  
#' @rdname registry
#' @export
packageRegistries <- function(regname=NULL, package=NULL, primary=FALSE){
	lns <- loadedNamespaces()
	if( !is.null(package) ) lns <- lns[lns %in% package]
	
	# early exit if no namespace
	if( !length(lns) ) return( character() )
	
	res <- lapply(lns, function(ns){
		reg <- packageRegistry(package=ns, quiet=TRUE)
		if( is.null(reg) ) return( character() )
		regnames <- reg$get_entry_names()
		res <- setNames(regnames, rep(ns, length(regnames)))
		if( primary ){
			pr <- sapply(res, function(n) reg$get_entry(n)$parent)
			res <- res[ nchar(pr) == 0L ]
		}
		res
	})

	res <- unlist(res)			
	if( !is.null(regname) ){
		res <- res[res == regname]
		if( primary && length(res) > 1L ){
			warning("Package registry - Found multiple primary registries '", regname, "' in packages "
					, str_out(res, Inf), " [using first one only]")
			res <- res[1L]
		}
	}
	res
}

#' \code{hasPackageRegistry} tells if a given package has a meta-registry or 
#' a given registry. 
#' 
#' @rdname registry
#' @export
hasPackageRegistry <- function(regname=NULL, package){
	isNamespaceLoaded(package) && !is.null( packageRegistry(regname, package=package, quiet=TRUE, entry=TRUE) )
}

#' @S3method format package_subregistry
format.package_subregistry <- function(x, ...){
	c(Key = x$key
	, Description = x$description
	, Entries = x$regobj$n_of_entries()
	, Parent = x$parent)	
}

#' @S3method format package_metaregistry
format.package_metaregistry <- function(x, ...){
	rec <- x$get_entries()
	data.frame(t(sapply(rec, base::format, ...))[, -1L, drop=FALSE])	 
}

#' @S3method print package_metaregistry
print.package_metaregistry <- function(x, ...){	
	NextMethod('print') #registry:::print.registry(x)
	print(format(x, ...))
}

#' @S3method xtable package_metaregistry 
#' @importFrom xtable xtable
xtable.package_metaregistry <- function(x, ...){
	d <- format(x)
	xtable::xtable(d, ...)
}

#' \code{setPackageRegistry} creates a package-specific registry within a package.
#'  
#' Each package sub-registry has its own set of fields.
#' Sub-registries defined by passing a character string in argument \code{regobj} of 
#' \code{setPackageRegistry} have the following fields: \code{'key'} and \code{'object'}
#' 
#' @param regobj a \code{\link[registry:regobj]{registry}} object or a single character 
#' string that indicates the class of the objects that are stored in the 
#' sub-registry.
#' See details for the list of the sub-registry's fields in this latter case.
#' @param description short description line about the registry.
#' It is recommended to provide such description as it makes clearer the purpose of the 
#' registry.
#' This description is shown when the registry object is printed/formated/listed.
#' @param entrydesc human readable description that is used in log messages 
#' when registering/removing entries.
#' @param ... named values used to set extra information about the new registry, that 
#' are stored in the corresponding fields of the meta-registry.
#' Currently not used, as no extra field other than \code{'description'} is defined.
#' @param overwrite a logical that indicate if an existing registry with the same 
#' should be overwritten if it exists.
#' 
#' @inheritParams packageRegistry
#' @rdname registry
#' @export
setPackageRegistry <- function(regname, regobj
								, description='', entrydesc=NA
								, ...
								, package=topenv(parent.frame())
								, overwrite=FALSE){
	
#	library(registry)
	
	# force overwrite in dev mode
	if( missing(overwrite) && isDevNamespace(package) ){
		overwrite <- TRUE
	}
	# check if sub-registry already exists
	oldreg <- packageRegistry(regname, quiet=TRUE, package=package)
	if( !is.null(oldreg) && !overwrite ){
		return( oldreg )
	}
	
	# get meta-registry (force creation)
	regenv <- .packageMetaRegistry(package, create=TRUE)
	nm <- packageSlot(regenv)
	ns_str <- str_c("package '", nm, "'")
	
	if( !is.null(oldreg) ){
		if( !overwrite ){
			if( isLoadingNamespace() ){ # exit if loading a namespace
				message("NOTE: Did not create registry '", regname,"' in ", ns_str, ": registry already exists.")
				return(oldreg)
			}
			stop("Could not create registry '", regname,"' in ", ns_str, ": registry already exists")
		}else{
			message("Remove registry '", regname,"' from ", ns_str)
			regenv$delete_entry(regname)
		}
	}
	message("Creating registry '", regname,"' in ", ns_str, ' ... ', appendLF=FALSE)
	
	.add_regclass <- function(x, newcl, before){
		cl <- class(x)
		ir <- which(cl == before)
		class(x) <- c(if( ir > 1 ) cl[1:(ir-1)] 
				, newcl, cl[ir:length(cl)])
		x
	}
	
	pkgregclass <- c(paste(regname, 'package_registry', sep='_'), 'package_registry')
	if( is.character(regobj) ){# regobj specifies the S4 class of the registry entries
		objtype <- regobj[1]
		regobj <- registry(entry_class = paste(regname, 'entry', sep='_')
						, registry_class = c(pkgregclass, 'object_subregistry'))
		# access key
		regobj$set_field("key", type="character", is_key = TRUE
				, index_FUN = match_partial_ignorecase)
		# object
		regobj$set_field("object", type=objtype, is_mandatory=TRUE, validity_FUN = validObject)
	}else if( is(regobj, 'registry') ){
		if( !is(regobj, 'package_registry') ){
			regobj <- .add_regclass(regobj, pkgregclass, 'registry')
		}
	}else{
		message('ERROR')
		stop("Invalid argument 'regobj': must be a class name or a registry object.")
	}
	# add field for REGISTERING package
	if( !"REGISTERINGpackage" %in% regobj$get_field_names() )
		regobj$set_field("REGISTERINGpackage", type='character', is_mandatory=TRUE, index_FUN=match_exact)
	# fix registry object
	regobj <- fix_registry(regobj)
	# add package
	attr(regobj, 'package') <- nm
	
	# create new meta entry
	regenv$set_entry(key=regname, regobj=regobj
					, description=description, entrydesc=entrydesc
					, ...)
	message('OK')
	# return newly created registry
	regenv$get_entry(regname)$regobj
}

regkeys <- function(regobj, ...){
	
	# get keys
	fld <- regobj$get_fields()
	keyfield <- names(fld[sapply(fld, function(x) isTRUE(x$is_key) )])
	if( nargs() == 1L ) return(keyfield)
	
	index_fields <- list(...)
	if( is.null(names(index_fields)) && length(index_fields)==1L )
		index_fields <- index_fields[[1L]]
	index_fields <- index_fields[!sapply(index_fields, is.null)]
	if( !length(index_fields) ) return(list())
	
	# remove fields that are not keys
	i <- match(keyfield, names(index_fields))
	index_fields[i[!is.na(i)]]
}

#' Finds an entry in a registry.
#' 
#' This function provides extra control on how entries are queried 
#' from a \code{\link[registry:regobj]{registry}} object.
#' 
#' @param regobj a registry object
#' @param ... key value(s) to look up.
#' If multiple indexes are used, then the primary key should come first.
#' @param all logical to indicate if hidden keys (starting with a '.') should be 
#' returned and output in message.
#' @param error a logical that indicates if an error should be thrown if the key has no match 
#' or multiple matches
#' @param exact a logical that indicates if matching should be exact or partial.
#' Note that if exact matches exist then they are returned, independently of the
#' value of \code{exact}.
#' @param KEYS alternative way of passing the key value(s).
#' If not missing, then arguments in \code{...} are discarded.
#' @param verbose a logical that indicates if verbosity should be toggle on
#' @param entry a logical that indicates if the 
#' @param msg a header to use in case of error.
#' 
#' @export
regfetch <- function(regobj, ..., all=FALSE, error=TRUE, exact=FALSE
						, KEYS = NULL
						, verbose=FALSE, entry=FALSE, msg=NULL){
	
	# load the registry package
#	library(registry)
	# list -- all -- keys if no key is specified
	keylist <- allkeys <- regobj$get_entry_names()
	if( !all ) keylist <- grep("^[^.]", keylist, value=TRUE)
	
	index_fields <- if( !is.null(KEYS) ){
		if( !is.list(KEYS) ) stop("Invalid argument <KEYS>: must be a list of field values.")
		KEYS
	}else list(...)
	# extract primary key
	key <- if( length(index_fields) ){
		# remove fields that are not keys if named list
		if( !is.null(names(index_fields)) )
			index_fields <- regkeys(regobj, index_fields)
		if( length(index_fields) ){
			paste(unlist(index_fields), collapse='_')
			str_out(index_fields, Inf, use.names=TRUE)
		}
	}
	if( is.null(key) ){
		return(keylist)
	}
	
	# set verbosity level
	if( !missing(verbose) ){
		ol <- lverbose(verbose)
		on.exit( lverbose(ol) )
	}
	
	if( !is.null(msg) ) msg <- str_c(msg, ' - ')
	
	if( regobj$n_of_entries() == 0L ){
		if( error )	stop(msg, "Registry is empty: no matching entry for key ", dQuote(key), ".")
		else return(NULL)
	}
	
	# get entry
	d <- do.call(regobj$get_entries, index_fields)
	# no entry found
	if( is.null(d) ){
		if( error ){
			stop(msg, "No matching entry for key ", dQuote(key), " in the registry."
							, "\n  Use one of: ", str_wrap(str_out(sort(allkeys), Inf), exdent=2), '.')
		}else return(NULL)
	}
	
	# look for exact matches
	if( is.list(index_fields) ){
		ex <- sapply(d, function(x) all(mapply(identical, index_fields, x[names(index_fields)])))
	}else{
		ex <- names(d) == index_fields
	} 
		
	# limit to exact mathes
	if( length(i <- which(ex)) ){
		d <- d[i]
	}else if( exact ){
		if( error ){
			stop(msg, "No exact match for key '", key, "' in the registry."
					, "\n  Use one of: ", str_wrap(str_out(allkeys, Inf), exdent=2), '.')
		}else return(NULL) 
	}
	
	if( all ) return(d)
	
	# multiple match
#	str(d)
	if( length(d) > 1L ){
		if( error ){
			stop(msg, "Multiple entries found for key ", dQuote(key), ": ", str_out(sort(names(d)), Inf), '.')
		}else return(NA)
	}
	
	# check single match
	if( length(d) != 1L )
		stop("Unexpected error: more than one entry was selected.")
	
	# return single match
	d <- d[[1L]]
	
	# return registry object if the entry is an automatic sub-registry
	if( !entry && is(regobj, 'object_subregistry') ) d$object
	else d
}
#' \code{pkgreg_fetch} fetches entries in a package registry, as set up by 
#' \code{\link{setPackageRegistry}}.
#' 
#' \code{pkgreg_fetch} loads the requested package registry and uses \code{regfetch} 
#' to retrieve data from it.
#' 
#' @inheritParams setPackageRegistry
#'  
#' @rdname regfetch
#' @export
pkgreg_fetch <- function(regname, ..., msg=NULL, where=topenv(parent.frame())){
	# get package registry
	regentry <- packageRegistry(regname, package=where, entry=TRUE, update=TRUE)
	# define addon error message
	if( missing(msg) && !is_NA(regentry$entrydesc) ) msg <- regentry$entrydesc
	# fetch from registry
	regfetch(regentry$regobj, ..., msg=msg)
}

#' \code{pkgreg_remove} removes an entry from a package registry.
#' 
#' @param quiet a logical that indicates if the operation should be performed quietly, 
#' without throwing errors or warnings.
#' 
#' @rdname regfetch
#' @export
pkgreg_remove <- function(regname, ..., msg=NULL, where=topenv(parent.frame()), quiet=FALSE){
	# get package registry
	regentry <- packageRegistry(regname, package=where, entry=TRUE, update=TRUE)
	# define addon error message
	if( missing(msg) && !is_NA(regentry$entrydesc) ) msg <- regentry$entrydesc
	# fetch from registry
	entry <- regfetch(regentry$regobj, ..., exact=TRUE, error=FALSE, all=TRUE, msg=msg)
	
	res <- if( !is.null(entry) ){
		# get the method registry and the method's fullname
		name <- names(entry)
		
		if( !quiet ){
			msg <- paste0("Removing ", msg, " '", name, "' from registry '", regname, "'")
			message(msg, ' ... ', appendLF=FALSE)
		}
		# delete from registry
		regentry$regobj$delete_entry(name)
		if( !quiet ) message('OK')
		TRUE
	}else{
		if( !quiet ){
			name <- str_out(list(...), Inf, use.names=TRUE)
			warning("Could not remove ", msg, " '", name, "': no matching registry entry.", call.=FALSE)
		}
		FALSE
	}

	if( quiet ) invisible(res)
	else res
}

extract_pkg <- function(x){
	sub("^(([^:]+)::)?(.*)", "\\2", x)
}

#' Automatic S4 Class for Registry Entries
#' 
#' @param registry a registry object
#' @param Class name of the class to generate
#' @param ... extra arguments passed to \code{\link{setClass}}.
#' 
setClassRegistry <- function(registry, Class, ...){
	
#	setClass(Class, representation, prototype, contains=character(),
#			validity, access, where, version, sealed, package,
#			S3methods = FALSE)
	
	f <- registry$get_fields()
	slots <- sapply(f, '[[', 'type', simplify=FALSE)
	
	args <- list(Class, representation=do.call('representation', slots))
#	if( !hasArg(validity) ){
#		.validity <-
#		sapply(f, function(x){
#			if(x$is_mandatory)
#				function(object){
#					if()
#				}
#		})
#		args$validity <- function(object){
#			
#		}
#	}
	do.call('setClass', c(args, ...))
}

#' \code{setPackageRegistryEntry} adds an entry in a package registry.
#' 
#' @param key entry identifier.
#' @param where package name or namespace that owns the registry. 
#' @param verbose a logical that indicates if verbosity should be toggle on.
#' @param msg addon message to print at the end of the output log line, 
#' when \code{verbose=TRUE}.
#' 
#' @rdname registry
#' @export
setPackageRegistryEntry <- function(regname, key, ..., overwrite=FALSE, verbose=FALSE
									, where=topenv(parent.frame()), msg=NULL){
	
	if( isLoadingNamespace() ){
		verbose <- TRUE
		if( missing(overwrite) ) overwrite <- TRUE
	}
	registry <- regname
	package <- where
	
	# check if the name provided is not empty
	if( nchar(key) == 0 ) stop('Invalid argument <key>: cannot be an empty string.')
	
	# build full key, that includes the name of the top calling namespace
	fullkey <- key
	top_ns <- topns(strict=FALSE)
	#
	
	# retrieve package registry (it must exist or this will throw an error)
	package <- packageEnv(package)
	subregentry <- packageRegistry(registry, package=package, entry=TRUE, update=TRUE)
	# get regobj (do that to ensure it is updated with entries from other packages)
	regobj <- subregentry$regobj
	
	# setup complete list of fields
	fields <- list(...)
	objdesc <- if( !is_NA(subregentry$entrydesc) ) subregentry$entrydesc else paste(registry, 'object')
	objdesc <- paste(objdesc, " '", key, "'", sep='')
	if( length(fields)==1L ){
		objdesc <- paste(objdesc, ' [', class(fields[[1L]]), ']', sep='')
		if( is.null(names(fields)) && is(regobj, 'object_subregistry') )
			names(fields) <- 'object'
	}
	fields$key <- key
	regpkg <- packageName(top_ns, .Global=TRUE)
	fields$REGISTERINGpackage <- regpkg
#	str(fields)
	#
	
	# check if the object is already registered
	oldentry <- regfetch(regobj, KEYS=fields, exact=TRUE, error=FALSE, all=TRUE)
	# error if already exists and not overwriting		
	if( !is.null(oldentry) && !overwrite ){ 
		if( verbose ) message("ERROR")
		stop("Cannot register ", objdesc, ": key already exists.")	
	}
	
	# add entry
	if( verbose ){
		action <- if( is.null(oldentry) ) 'Registering' else 'Replacing'
		message(action, " ", objdesc, msg, " ... ", appendLF=FALSE)
	}
	# delete old entry
	if( !is.null(oldentry) ){
		regobj$delete_entry(names(oldentry)[1L])
	}
	# do add entry
	do.call(regobj$set_entry, fields)
	if( verbose ) message("OK")
	
	# if the registration happens during loading another package: 
	# create local registry and add entry to it.
	# It will be merged to the main registry on the next call to 
	# packageRegistry after the package is loaded.
	lns <- getLoadingNamespace(env=TRUE)
	if( !is.null(lns <- getLoadingNamespace(env=TRUE)) && !identical(lns, package) ){
		# clone registry
		if( nchar(subregentry$parent) ){
			warning("Deriving package registry '", registry, "' in package ", lns
					, " from ", subregentry$parent, " instead of ", subregentry$package, immediate.=TRUE)
			parent <- subregentry$parent
		}else parent <- subregentry$package
		fullregistry <- str_c(parent, '::', registry)
		
		if( is.null(locregobj <- packageRegistry(fullregistry, package=lns, quiet=TRUE)) ){
			# clone registry
			locregobj <- clone_regobj(regobj, empty=TRUE)
			# attach to loading namespace
			locregobj <- setPackageRegistry(fullregistry, locregobj
											, description = subregentry$description
											, entrydesc = subregentry$entrydesc
											, parent = parent
											, package = lns)
		}
		
		action <- 'Adding'
		if( !is.null(locentry <- regfetch(locregobj, KEYS=fields, exact=TRUE, error=FALSE, all=TRUE)) ){
			action <- 'Overwriting'
			locregobj$delete_entry(names(locentry)[1L])
		}
		# add entry into local registry
		if( verbose ) message(action, " entry '", key, "' in registry '", fullregistry, "' ... ", appendLF=FALSE)
		do.call(locregobj$set_entry, fields)
		if( verbose ) message("OK")
	}
	#
	
	# return registered object
	regfetch(regobj, KEYS=fields, exact=TRUE)

}

# clone a registry object
clone_regobj <- function(regobj, empty=FALSE){
	tmp <- tempfile('registry')
	on.exit(unlink(tmp))
	saveRDS(regobj, file=tmp)
	newreg <- readRDS(tmp)
	# empty entries if necessary
	if( empty ){
		sapply(newreg$get_entry_names(), newreg$delete_entry)
	}
	newreg
}
