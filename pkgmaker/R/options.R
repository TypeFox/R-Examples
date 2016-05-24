# Package specific option system
# 
# Author: Renaud Gaujoux
# Creation: 25 Apr 2012
###############################################################################

#' @include unitTests.R
#' @include devutils.R
NULL

#' Package Specific Options
#' 
#' The following functions to access/set the options from the set are assigned 
#' in \code{envir}:
#' \describe{
#' \item{<subset>Options}{}
#' \item{<subset>GetOption}{}
#' }
#' 
#' @param ... a single named list or named arguments that provide the default 
#' options and their values.
#' @param NAME name of the set of options.
#' This is used as a prefix for the name of the associated global 
#' option: \code{package:<name>}.
#' @param ENVIR environment where the option wrapper functions will be defined.
#' No function is defined if \code{ENVIR=NULL} 
#' @param RESET a logical that indicates whether the option set should overwrite
#' one that already exists if necessary. 
#' The default is \code{FALSE} (i.e. no reset), except when loading a namespace, 
#' either from an installed package or a development package -- with devtools. 
#' If \code{FALSE}, an error is thrown if trying to setup options with the same name.
#'
#' @export
setupPackageOptions <- function(..., NAME=NULL, ENVIR=topenv(parent.frame()), RESET = isLoadingNamespace()){
	
	defaults <- .list_or_named_dots(...)
	
	# do not write into the Global environment
	e <- parent.frame()
	if( missing(ENVIR) && identical(e, .GlobalEnv) ){
		ENVIR <- NULL
	}
	
	# get calling package
	pkg <- packageName(.Global=TRUE)
	
	# prefix for the wrapper functions
	fprefix <- if( is.null(NAME) ) tolower(pkg) else NAME
	
	# define name for the option set
	optname <- pkg
	if( !is.null(NAME) )
		optname <- paste(optname, NAME, sep=':')
	
	# create package_options object
	optobj <- as.package_options(optname, defaults=defaults)
	
	# check if options with the same key are not already registered
	OLD <- getOption(optobj$name)
	if( !is.null(OLD) && !RESET )
		stop("Package specific options '", OLD$name, "' already exist: " 
				, " (", length(OLD$options())," default option(s))")
	
	# register the package_options object in global options
	message(if( is.null(OLD) ) "Setting" else "Resetting"
			, " package specific options: ", optobj$name
			, " (", length(optobj$options())," default option(s))")
	options(setNames(list(optobj), optobj$name))
	
	# (re)load registered package_options object from global options
	optobj <- getOption(optobj$name)
	stopifnot( !is.null(optobj) )
	
	# define wrapper functions in the supplied environment
	if( !is.null(ENVIR) ){
		isfun <- unlist(eapply(optobj, is.function))
		isfun <- isfun[names(isfun) != 'newOptions']
		ifun <- which(isfun)
		lapply(names(isfun)[ifun], function(x){
			f <- get(x, envir=optobj)
			assign(paste(fprefix, x, sep='.'), f, envir=ENVIR)
		})
	}
	
	# return package_options object
	optobj
}

is.package_options <- function(x){
	is(x, 'package_options')
}

#' @S3method print package_options
print.package_options <- function(x, ...){
	cat("<Package specific options: ", x$name, ">\n", sep='')
	cat("Registered: ", !is.null(getOption(x$name)), "\n", sep='')
	def <- if( identical(x$.options, x$.defaults) ) " <as default>"
	# show options
	if( length(x$.options) ){
		cat("Options",def,":\n", sep='');
		str(x$.options) 
	}else 
		cat("Options: none\n")
	# show defaults
	if( is.null(def) ){
		cat("Defaults:\n"); str(x$.defaults)
	}
}


#' \code{option_symlink} creates a symbolic link to option \code{x}.
#' 
#' @export
#' @rdname options
option_symlink <- function(x){
	if( !is.character(x) )
		stop("Symbolic link options must be character strings")
	structure(x, class='option_symlink')
}
#' \code{is_option_symlink} tests if \code{x} is a symbolic link option.
#' 
#' @param opts a list of options
#' 
#' @export
#' @rdname options
is_option_symlink <- function(x, opts){
	if( missing(opts) ) is(x, 'option_symlink')
	else is(opts[[x]], 'option_symlink')
}

#' \code{option_symlink_target} returns the end target option of a symbolic link 
#' option \code{x}.
#' 
#' @export
#' @rdname options
option_symlink_target <- function(x, opts){
	
	if( !is.list(opts) )
		stop("invalid argument `opts`: must be a list object")
	
	n <- 0
	track <- NULL
	while( is_option_symlink(x, opts) ){
		if( x %in% track )
			stop("cycling symbolic link options: ", str_out(c(track, x), Inf, sep=' -> '))
		track <- c(track, x)
		x <- opts[[x]]
		n <- n + 1
		
	}
	x
	
}

# unit test for option symbolic links
unit.test('option_symlink', {

	opt <- setupPackageOptions(a=1,b=2,c=option_symlink('a'),d=4)
	
	.test <- function(msg){
		checkIdentical(names(opt$options('a')), 'a', paste(msg, " - options: name of target is ok"))
		checkIdentical(names(opt$options('c')), 'c', paste(msg, " - options: name of link is ok"))
		checkIdentical(opt$options('c'), setNames(opt$options('a'), 'c'), paste(msg, " - options: link ok"))
		checkIdentical(opt$getOption('a'), opt$getOption('c'), paste(msg, " - getOption: link ok"))
	}
	
	.test('Default')
	opt$options(a=100)
	.test('After setting target')
	opt$options(c=50)
	.test('After setting link')
			
})

#' \code{as.package_options} creates an object such as the 
#' ones used to stores package specific options.
#' 
#' @param x a character string, a list or an object of class 
#' \code{package_options}.
#' @param defaults \code{NULL} or a list of default options 
#' with their values.   
#'
#' @export
#' @rdname options
as.package_options <- function(..., defaults=NULL){
	
	args <- .list_or_named_dots(...)
	
	x <- if( is.null(names(args)) ) args[[1]] 
	if( !is.null(names(args)) ) defaults <- args
	if( is.null(x) ) x <- basename(tempfile(''))
	
	# early exit if already a package_options object
	if( is.package_options(x) ){
		
		# new defaults?: clone into a new package_options object
		if( !missing(defaults) && is.list(defaults) ){
			optname <- basename(tempfile(str_c(x$name, '_')))
			x <- as.package_options(x$.options, defaults)
			x$name <- optname
		}
	
		return(x)
	}
	
	# create a package_options object
	.OPTOBJ <- structure(list2env(list(name=NULL, .options=NULL, .defaults=defaults))
						, class='package_options')
	
	if( is.character(x) ){
		
		# build name as 'package:*'
		x <- sub("^package:", '', x)
		.OPTOBJ$name <- paste('package:', x[1L], sep='')
		
	}else if( is.list(x) ){
		.OPTOBJ$name <- tempfile('package:')
		.OPTOBJ$.options <- x
	}else
		stop("Invalid argument `x`: must be a character string or a list.")
	
	# define options() 
	.OPTOBJ$options <- function(...){
		# call .options on package_options object
		.options(..., .DATA=.OPTOBJ)
	}
	# define getOption
	.OPTOBJ$getOption <- function (x, default = NULL) 
	{
		# use local specific function options()
		options <- .OPTOBJ$options
		
		if (missing(default)) 
			return(options(x)[[1L]])
		if (x %in% names(options())) 
			options(x)[[1L]]
		else default
	}
	# define newOption
	.OPTOBJ$newOptions <- function(...){
		defs <- .list_or_named_dots(..., named.only=TRUE)
		
		lapply(seq_along(defs),
			function(i){
			name <- names(defs)[i]
			value <- defs[[i]]
			# check defaults
            in_opts <- name %in% names(.OPTOBJ$.defaults) && !identical(.OPTOBJ$.defaults[[name]], value)
			if( in_opts && !isLoadingNamespace() ){
				message("Skipping option ", .OPTOBJ$name, "::`", name, "`: already defined with another default value")
		    }else{
                if( in_opts )
                    message("Overwriting option ", .OPTOBJ$name, "::`", name, "` : already defined with another default value")
				.OPTOBJ$.defaults[[name]] <- value
				.OPTOBJ$.options[[name]] <- value
			}
		})
		invisible()
	}
	# define resetOptions
	.OPTOBJ$resetOptions <- function(..., ALL=FALSE){
		
		defaults <- .OPTOBJ$.defaults
		if( ALL ){
			.OPTOBJ$.options <- NULL
		}
		if( length(list(...)) > 0L ){
			onames <- c(...)
			if( !is.character(onames) )
				stop('character strings expected for resetting option names')
			defaults <- defaults[names(defaults) %in% onames]
			if( length(not_default <- onames[!onames %in% names(defaults)]) ){
				.OPTOBJ$.options[not_default] <- NULL
			}
		}
		if( length(defaults) ){
			.OPTOBJ$options(defaults)
		}
	}
	# define showOptions
	.OPTOBJ$printOptions <- function() print(.OPTOBJ)
	
	# initialise with default options 
	.OPTOBJ$resetOptions()
	
	# return pacakge_options object
	.OPTOBJ
}


#' The method \code{[[} is equivalent to \code{options()} or \code{getOption(...)}:
#' e.g. \code{obj[[]]} returns the list of options defined in \code{obj}, and 
#' \code{obj[['abc']]} returns the value of option \code{'abc'}.
#' 
#' @param ... arguments passed to \code{getOption} (only first one is used). 
#'  
#' @rdname options 
#' @S3method [[ package_options
"[[.package_options" <- function(x, ...){
	if( missing(..1) ) x$options()
	else x$getOption(..1)
}

#' @S3method [[<- package_options
"[[<-.package_options" <- function(x, i, value){
	x$.options[[i]] <- value 
}


##' @S3method [[ package_options
#`[[.package_options` <- function(x, ..., follow=FALSE){
#	
#	if( missing(..1) ) as.list(x$.options)
#	else if( follow ){
#		x$.options[[option_symlink_target(..1, x)]]
#	}else x$.options[[..1]]
#}
#
##' @S3method [[<- package_options
#`[[<-.package_options` <- function(x, i, ..., value){
#	
#	follow <- if( missing(..1) ) FALSE else ..1 
#	if( follow ){
#		old <- x[[i]]
#		if( is_option_symlink(old) && !is_option_symlink(value) )
#			x$.options[[option_symlink_target(i, x)]] <- value
#	}else x$.options[[i]] <- value
#}

.list_or_named_dots <- function(..., named.only=FALSE){
	
	dots <- list(...)
	if( length(dots) == 0L ) return()
	
	params <- dots
	if( is.null(names(dots)) && length(dots)==1L ){
		if ( is.list(dots[[1L]]) ){ 
			params <- dots[[1L]]
			if( is.null(names(params)) || any(names(params)=='') )
				stop("single list argument must only have named elements")
		}
	}
	if( named.only ){
		if( is.null(names(params)) || any(names(params)=='') )
			stop("all arguments be named")
	}
	
	params
}

#' Quick Option-like Feature
#' 
#' \code{mkoptions} is a function that returns a function that 
#' behaves like \code{\link[base]{options}}, with an attached 
#' internal/local list of key-value pairs.
#' 
#' @rdname local-options
#' @seealso \code{\link{setupPackageOptions}}
#' @export
#' 
#' @examples
#' f <- mkoptions(a=3, b=list(1,2,3))
#' str(f())
#' f('a')
#' f('b')
#' str(old <- f(a = 10))
#' str(f())
#' f(old)
#' str(f())
#' 
mkoptions <- function(...){
	
	.DATA <- new.env(parent=emptyenv())
	.defaults <- list(...)
	.DATA$.options <- list(...)
	function(...){		
		.options(..., .DATA=.DATA)
	}
}

#' \code{.options} is a low-level function that mimics the behaviour 
#' of the base function \code{\link[base]{options}}, given a set 
#' of key-value pairs.
#' It is the workhorse function used in \code{mkoptions} and package-specific
#' option sets (see \code{\link{setupPackageOptions}})
#' 
#' @param ... list of keys or key-value pairs.
#' For \code{mkoptions} these define inital/default key-value pairs. 
#' @param .DATA a list or an environment with an element \code{.options}.
#' 
#' @rdname local-options
.options <- function(..., .DATA){
	
	opts <- if( is.package_options(.DATA) || is.environment(.DATA) ) .DATA$.options else .DATA
	
	params <- .list_or_named_dots(...)
	# return complete option list if no other argument was passed
	if( is.null(params) ) return(opts)
	
	# initialise opts to an empty list if necessary 
	if( is.null(opts) ) opts <- list()
	stopifnot( is.list(opts) )
	
	# READ ACCESS
	if ( is.null(names(params)) ){
		if( !is.character(c(...)) )
			stop('character strings expected for option names')
		
		cparams <- c(...)
		# retrieve options as a list (use sapply to also get non-existing options)
		res <- sapply(cparams, function(n){
					# follow link if necessary
					opts[[option_symlink_target(n, opts)]]
				}, simplify=FALSE)
		return(res)
	}
	
	# WRITE ACCESS
	old <- sapply(names(params), 
			function(name){
				# assign the new value into the options environment
				val <- params[[name]]
				old <- opts[[name]]
				# change value of target if symlink and the new value is not a symlink
				if( is_option_symlink(old) && !is_option_symlink(val) )
					opts[[option_symlink_target(name, opts)]] <<- val
				else
					opts[[name]] <<- val
				# return the option's old value
				old
			}
			, simplify = FALSE
	)	
	#old <- old[!sapply(old, is.null)]

	# update package_options object in place if necessary (NB: it is an environment)
	if( is.package_options(.DATA) || is.environment(.DATA) ) .DATA$.options <- opts
	
	# return old values of the modified options
	return( invisible(old) )
}

#' \code{packageOptions} provides access to package specific options from a 
#' given package that were defined with \code{setupPackageOptions}, and behaves as the base function \code{\link[base]{options}}.
#' 
#' @param PACKAGE a package name
#' @inheritParams base::options 
#' 
#' @export
#' @rdname options
packageOptions <- function(..., PACKAGE = packageName()){
		
	# create/retrieve a package_options object from .DATA
	optobj <- as.package_options(PACKAGE)
	optobj <- getOption(optobj$name)
	
	# call the package_options object's options() function
	optobj$options(...)
}

#' \code{listPackageOptions} returns the names of all option 
#' currently defined with \code{setupPackageOptions}.
#' 
#' @return a character vector (possibly empty).
#'
#' @export
#' @rdname options 
#' @examples
#' listPackageOptions()
#' 
listPackageOptions <- function(){
	grep('^package:', names(options()), value=TRUE)
} 
