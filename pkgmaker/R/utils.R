# General utility functions
# 
# Author: Renaud Gaujoux
# Creation: 25 Apr 2012
###############################################################################

# or-NULL operator (borrowed from Hadley Wickham)
'%||%' <- function(x, y) if( !is.null(x) ) x else y

#' Get Anywhere
#' 
#' Similar to \code{\link{getAnywhere}}, but looks for the value of its argument. 
#' 
#' @param x a single character string
#' 
#' @export
cgetAnywhere <- function(x){
	do.call("getAnywhere", list(x))
}

#' Silent Require
#' 
#' Silently require a package.
#' 
#' @inheritParams base::require
#' @param ... extra arguments passed to \code{\link{require}}.
#' 
#' @export
require.quiet <- function(package, character.only = FALSE, ...){
	
	if( !character.only )
		package <- as.character(substitute(package))
	utils::capture.output(suppressMessages(suppressWarnings(
	 res <- do.call('require', 
			 list(package=package, ..., character.only=TRUE, quietly=TRUE))
	)))
	res
}


#' Testing R Version
#' 
#' Compares current R version with a given target version, which may be useful  
#' for implementing version dependent code. 
#' 
#' @param x target version to compare with.
#' @param test numeric value that indicates the comparison to be carried out.
#' The comparison is based on the result from 
#' \code{utils::compareVersion(R.version, x)}:
#' \itemize{
#' \item 1: is R.version > \code{x}?
#' \item 0: is R.version = \code{x}?
#' \item -1: is R.version < \code{x}?
#' } 
#' 
#' @return a logical
#' @export
#' @examples
#' 
#' testRversion("2.14")
#' testRversion("2.15")
#' testRversion("10")
#' testRversion("10", test = -1)
#' testRversion("< 10")
#' testRversion(Rversion())
#' testRversion(paste0('=', Rversion()))
#' 
testRversion <- function(x, test=1L){
	rv <- Rversion()
    op <- '=='
    if( grepl("^[=<>]", str_trim(x)) ){
        m <- str_match(x, "^([<>=]=?)(.*)")
        if( is.na(m[, 1]) ) stop('Invalid version specification: ', x)
        op <- m[, 2]
        if( op == '=' ) op <- '=='
        x <- str_trim(m[, 3L])
        if( !missing(test) ) warning("Ignoring argument `test`: comparison operator was passed in argument `x`")
        test <- 0L
    }
	do.call(op, list(utils::compareVersion(rv, x), test))
}

#' Complete R version
#' 
#' Returns the complete R version, e.g. 2.15.0
#' 
#' @export
#' @examples
#' Rversion()
#' 
Rversion <- function(){
	paste(R.version$major, R.version$minor, sep='.')
}

#' Formatting Utilities
#' 
#' \code{str_out} formats character vectors for use in show methods or 
#' error/warning messages.
#' 
#' @param x character vector
#' @param max maximum number of values to appear in the list. If \code{x} has 
#' more elements than \code{max}, a \code{"..."} suffix is appended.
#' @param quote a logical indicating whether the values should be quoted with 
#' single quotes (defaults) or not. 
#' @param use.names a logical indicating whether names should be added to the 
#' list as \code{NAME=VAL, ...} or not (default).
#' @param sep separator character
#' @param total logical that indicates if the total number of elements should be 
#' appended to the formatted string as \code{"'a', ..., 'z' (<N> total)"}.  
#' 
#' @return a single character string
#' 
#' @examples
#' 
#' x <- letters[1:10]
#' str_out(x)
#' str_out(x, 8)
#' str_out(x, Inf)
#' str_out(x, quote=FALSE)
#' str_out(x, total = TRUE)
#' 
#' @export
str_out <- function(x, max=3L, quote=is.character(x), use.names=FALSE, sep=", ", total = FALSE){
	if( is_NA(max) ) max <- Inf
	suffix <- NULL
    nTotal <- length(x)
	if( max > 2 && length(x) > max ){
		suffix <- "..."
		x <- c(head(x, max-1), tail(x, 1))
	}
	x <- head(x, max)
	
	# add quotes if necessary
	quote <- 
			if( isTRUE(quote) ) "'"
			else if( is.character(quote) ) quote
	if( !is.null(quote) ) x <- unlist(lapply(x, function(v) paste(quote,v,quote, sep='')))
	else if( all(sapply(x, isInteger)) ) x <- unlist(lapply(x, function(v) str_c(v,'L')))
	# add names if necessary
	if( use.names && !is.null(names(x)) ){
		nm <- str_c(names(x),'=')
		x <- paste(ifelse(nm=='=','',nm), x, sep='')
	}
	# insert suffix
	if( !is.null(suffix) ){
		x <- c(head(x, length(x)-1L), suffix, tail(x, 1L))
	}
	s <- paste(paste(x, collapse=sep), sep='')
	
	if( total ) s <- paste0(s, ' (', nTotal, ' total)')
	
	# return formatted string 
	s
}

#' \code{str_desc} builds formatted string from a list of complex values.
#' 
#' @param object an R object
#' @param exdent extra indentation passed to str_wrap, and used if the output 
#' should spread over more than one lines.
#' 
#' @rdname str_out
#' @export
str_desc <- function(object, exdent=0L){
	p <- sapply(object, function(x){
				if( is.atomic(x) && length(x) == 1L ) x
				else paste("<", class(x), ">", sep='')
			})
	str_wrap(str_out(p, NA, use.names=TRUE, quote=FALSE), exdent=exdent)
}

#' \code{str_fun} extracts and formats a function signature.
#' It typically formats the output \code{capture.output(args(object))}.
#' @rdname str_out
#' @export
#' @examples 
#' str_fun(install.packages)
str_fun <- function(object){
	s <- capture.output(args(object))
	paste(s[-length(s)], collapse="\n")
}

# From example in ?toupper
capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s,1,1)),
                {s <- substring(s,2); if(strict) tolower(s) else s},
                sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

#' Finding Differences Between Strings
#' 
#' Computes which characters differ between two strings.
#' 
#' @param x a single string
#' @param y a single string
#' @return an integer vector containing the index of all mis-matched characters
#' in the first string.
#' @export
#' 
#' @examples
#' 
#' # strings to compare
#' x <- "once upon a time"
#' y <- "once upon a time there was"
#' z <- "once upon two times"
#' 
#' # diff: x - y
#' d <- str_diff(x, y)
#' d
#' str(d)
#' 
#' # other comparisons 
#' str_diff(y, x)
#' str_diff(x, x)
#' str_diff(x, z)
#' str_diff(y, z)
#' 
str_diff <- function(x, y){
	sx <- strsplit(x,'')[[1]]
	sy <- strsplit(y,'')[[1]]
	n <- min(length(sx), length(sy))
	res <- mapply('!=', head(sx,n), head(sy,n))
	wres <- which(res)
	if( length(sx) > length(sy) )
		wres <- c(wres, (n+1):length(sx))
	attr(wres, 'str') <- list(x=x,y=y)
	class(wres) <- 'str_diff'
	wres
}

#' @S3method print str_diff
print.str_diff <- function(x, ...){
	s <- attr(x, 'str')
	n <- max(nchar(s$x), nchar(s$y))
	d <- rep('.', n)
	d[x] <- '*'
	if( (n2 <- nchar(s$y)-nchar(s$x)) )
		d[(n-abs(n2)+1):n] <- if( n2 > 0L ) '-' else '+'
	cat(str_c(s$x, collapse=''), "\n")
	cat(str_c(d, collapse=''), "\n")
	cat(str_c(s$y, collapse=''), "\n")				
}

#' Extracting Local Function Definition
#' 
#' @description
#' \code{extractLocalFun} Extracts local function from wrapper functions of the following type, typically 
#' used in S4 methods:
#' \samp{
#' function(a, b, ...)\{
#' 	.local <- function(a, b, c, d, ...)\{\}
#'	.local(a, b, ...)
#' \}
#' }
#'
#' @param f definition of the wrapper function
#' 
#' @return a function
#' @export
#' @rdname formals
extractLocalFun <- function(f){
	bf <- body(f)
	
	txt <- as.character(bf)[2]
	# in R-2.14.2 -- at least, as.character does not return the complete body
	# so some text manipulation is necessary 
	if( !grepl("\\{", txt) ){
		sf <- capture.output(print(bf))
		w <- tail(grep("^\\s*\\.local\\(", sf), 1L)
		txt <- paste(sf[-w], collapse="\n")
	}
	expr <- parse(text=txt)
	e <- new.env()
	eval(expr, e)
} 

#' Extended Formal Extraction
#'
#' Works for methods that are created (setMethod) as a wrapper function to an 
#' internal function named .local.
#'
#' @inheritParams extractLocalFun
#' @return a paired list like the one returned by \code{\link{formals}}. 
#' 
#' @export
#' @import codetools
#' @rdname formals
allFormals <- function(f){
	
	# look inside method for S4 methods
	if( is(f, 'MethodDefinition') ){
		
		# check if the method is defined as a wrapper function
		f <- f@.Data
		lf <- try(codetools::getAssignedVar(body(f)), silent=TRUE)
		if( !identical(lf, '.local') ) return( formals(f) )
		# extract arguments from local function
		lfun <- extractLocalFun(f)
		res <- formals(lfun)
		# set default values from the generic, only for arguments that have no 
		# default values in the method
		generic_args <- formals(f)
		meth_no_default <- sapply(res, is.symbol) 
		gen_no_default <- sapply(generic_args, is.symbol)
		generic_args <- generic_args[ !gen_no_default ]
		generic_args <- generic_args[ names(generic_args) %in% names(res[meth_no_default]) ]
		if( length(generic_args) ){
			res[names(generic_args)] <- generic_args
		}
		# return complete list of arguments
		res
		
	}else if( is.function(f) ) formals(f)
	
}

#' Alternative S4 Constructor
#' 
#' An alternative version of \code{\link{new}} to create objects based on a list
#' of values. 
#' 
#' @param class Class name to instanciate
#' @param ... extra arguments from which slot values are extracted by exact 
#' matching of names.
#' 
#' @export
#' @examples
#' 
#' setClass('A', contain='character', representation(x='numeric', y='character'))
#' 
#' # identical behaviour with standard calls
#' identical(new('A'), new2('A'))
#' identical(new('A', x=1), new2('A', x=1))
#' 
#' # but if passing that are names not slots 
#' identical(new('A'), new2('A', b=1))
#' identical(new('A', x=1), new2('A', x=1, b=3))
#' identical(new('A', x=1), new2('A', x=1, b=3))
#' 
#' # standard `new` would coerce first unnamed argument into parent of 'A' (i.e. 'character') 
#' new('A', list(x=1))
#' new('A', list(x=1, y='other'))
#' # `new2` rather use it to initialise the slots it can find in the list 
#' identical(new('A', x=1), new2('A', list(x=1)))
#' identical(new('A', x=1, y='other'), new2('A', list(x=1, y='other')))
#' 
#' 
new2 <- function(class, ...){
	sl <- getSlots(class)
	if( nargs() == 1L ) return( new(class) )
	
	dots <- list(...)
	if( nargs() == 2L && is.null(names(dots)) ){
		l <- dots[[1]]
		if( !is.list(l) )
			stop("Invalid call: single unnamed argument must be a list")
		dots <- l
	}
	
	if( is.null(names(dots)) || any(names(dots)=='') )
		stop("Invalid call: all slot arguments must be named")
	dots <- dots[names(dots) %in% names(sl)]
	do.call('new', c(list(class), dots))
}


#' One-off Global Variables
#' 
#' Defines a function that allow to get/assign a global variable whose value is 
#' ensured to be reset after each access.
#'   
#' @param default default value to which the global variable is reset after each 
#' access. Default is \code{NULL}.
#' 
#' @return a function with one argument (\code{value}) that provides get/set access
#' to a global variable.
#' If called with a value, it assigns this value to the global variable.
#' If called with no argument, it returns the current value of the global variable and 
#' reset it to its default value -- as defined at its creation. 
#'
#' @export
#' 
#' @examples
#' 
#' x <- oneoffVariable(0)
#' # returns default value
#' x()
#' # assign a value
#' x(3)
#' # get the value
#' x()
#' # second call returns default value again 
#' x()
#'  
oneoffVariable <- function(default=NULL){
	.var <- default
	function(value){
		if( missing(value) ){
			res <- .var
			.var <<- default
			res
		}else
			.var <<- value
	}
}


##  Exit Error Checker
##  
##  This function defines a function that checks if an error has been 
##  thrown after its definition.
##  It may be used to perform tasks on function exit depending on 
##  how the function exit (normal return or with an error).
##  
##  The function \code{errorCheck} itself is meant to be called at 
##  the beginning of functions that use \code{\link{on.exit}} to 
##  perform tasks when exiting.
##  The error checker function returned, when used in \code{on.exit} 
##  expressions, enables to distinguish between a normal exit and 
##  an exit due to an error, allowing is to perform tasks specific 
##  to each scenario.
##  
##  IMPORTANT: this function is not 100\% perfect in the sense that 
##  it will detect an error as soon as one has been thrown, even it 
##  is catched before the exit -- with \code{\link{try}} or 
##  \code{\link{tryCatch}}.
##  
##  @export
##  @examples 
##  
##  # define some function
##  f <- function(err){
##  
##   # initialise an error checker
##  	isError <- errorCheck()
##  
##   # do something on exit that depends on the error status
##  	on.exit({
##  		if(isError()) cat("with error: cleanup\n") 
##  		else cat("no error: do nothing\n") 
##  	})
##  	
##   # throw an error here
##  	if( err ) stop('There is an error')
##   
##  	1+1
##  }
##  
##  # without error
##  f(FALSE)
##  # with error
##  try( f(TRUE) )
##  
#errorCheck <- function(){
#	
#	# initialise with unique error message
#	.err <- tryCatch(stop('ERROR_CHECK:', digest(tempfile())), error=function(e) conditionMessage(e))
#	tb_digest <- function() digest(capture.output(traceback(max.lines=NULL)))
#	.traceback <- tb_digest()
#	
#	function(){
#		# error message is different
#		# tb_digest() != .traceback
#		length(grep(.err, msg, fixed=TRUE, invert=TRUE)) == 1L
#	}
#}


#' Global Static Variable
#' 
#' \code{sVariable} defines a function that acts as a global
#' static variable.
#' 
#' @param default default value for the static variable. 
#' 
#' @export
#' @examples 
#' 
#' # define variable
#' x <- sVariable(1)
#' # get value (default)
#' x()
#' # set new value: return old value
#' old <- x(3)
#' old
#' # get new value
#' x()
#' 
sVariable <- function(default=NULL){
	.val <- default
	function(value){
		if( missing(value) ) .val
		else{
			old <- .val
			.val <<- value
			old
		}
	}
}

#' Exit Error Checks
#' 
#' \code{exitCheck} provides a mechanism to distinguish the exit status
#' in \code{\link{on.exit}} expressions.
#' 
#' It generates a function that is used wihtin a function's body to 
#' "flag" normal exits and in its \code{\link{on.exit}} expression
#' to check the exit status of a function.
#' Note that it will correctly detect errors only if all normal exit 
#' are wrapped into a call to it. 
#' 
#' @export
#' 
#' @examples
#' 
#' # define some function
#' f <- function(err){
#' 
#'  # initialise an error checker
#' 	success <- exitCheck()
#' 
#'  # do something on exit that depends on the error status
#' 	on.exit({
#' 		if(success()) cat("Exit with no error: do nothing\n") 
#' 		else cat("Exit with error: cleaning up the mess ...\n") 
#' 	})
#' 	
#'  # throw an error here
#' 	if( err ) stop('There is an error')
#'  
#' 	success(1+1)
#' }
#' 
#' # without error
#' f(FALSE)
#' # with error
#' try( f(TRUE) )
#' 
exitCheck <- function(){
	
	.success <- FALSE
	function(x){
		if( nargs() == 0L ) .success
		else{
			.success <<- TRUE
			x
		}
	}
}

#' Ordering Version Numbers
#' 
#' Orders a vector of version numbers, in natural order.
#' 
#' @param x a character vector of version numbers
#' @param decreasing a logical that indicates if the ordering should be decreasing
#' 
#' @export
#' @examples
#' 
#' v <- c('1.0', '1.03', '1.2')
#' order(v)
#' orderVersion(v)
#' 
orderVersion <- function(x, decreasing=FALSE){
	tx <- gsub("[^0-9]+",".", paste('_', x, sep=''))
	stx <- strsplit(tx, ".", fixed=TRUE)
	mtx <- max(sapply(stx, length))
	tx <- sapply(stx, 
			function(v) paste(sprintf("%06i", c(as.integer(v[-1]),rep(0, mtx-length(v)+1))), collapse='.')
	)	
	order(tx, decreasing=decreasing)
}

#' @param ... extra parameters passed to \code{orderVersion}
#' 
#' @export
#' @rdname orderVersion
#' @examples
#' 
#' sort(v)
#' sortVersion(v)
sortVersion <- function(x, ...){
	x[orderVersion(x, ...)]
}

#' Checking for Missing Arguments
#' 
#' This function is identical to \code{\link{hasArg}}, except that 
#' it accepts the argument name as a character string.
#' This avoids to have a check NOTE about invisible binding variable.  
#' 
#' @param name the name of an argument as a character string.
#' 
#' @export
#' @examples
#' 
#' f <- function(...){ hasArg2('abc') }
#' f(a=1)
#' f(abc=1)
#' f(b=1)
#' 
hasArg2 <- function (name) 
{
	name <- as.name(name)
	## apply methods::hasArg
	aname <- as.character(substitute(name))
	fnames <- names(formals(sys.function(sys.parent())))
	if (is.na(match(aname, fnames))) {
		if (is.na(match("...", fnames))) 
			FALSE
		else {
			dotsCall <- eval(quote(substitute(list(...))), sys.parent())
			!is.na(match(aname, names(dotsCall)))
		}
	}
	else eval(substitute(!missing(name)), sys.frame(sys.parent()))
	##
}

#' Exposing Object Attributes
#' 
#' The function \code{ExposeAttribute} creates an S3 object that 
#' exposes all attributes of any R object, by making them accessible via 
#' methods \code{\link{$}} and/or \code{\link{$<-}}.
#' 
#' @param object any R object whose attributes need to be exposed
#' @param ... attributes, and optionally their respective values or
#' access permissions.
#' See argument \code{value} of \code{attr_mode} for details on the
#' way of specifying these. 
#' @param .MODE access mode:
#' \describe{
#' \item{\dQuote{r}:}{ (read-only) only method \code{$} is defined}
#' \item{\dQuote{w}:}{ (write-only) only method \code{$<-} is defined}
#' \item{\dQuote{rw}:}{ (read-write) both methods \code{$} and \code{$<-} 
#' are defined}
#' } 
#' @param .VALUE logical that indicates if the values of named arguments 
#' in \code{...} should be considered as attribute assignments, 
#' i.e. that the result object has these attributes set with the specified values.
#' In this case all these attributes will have the access permission
#' as defined by argument \code{.MODE}.
#' 
#' @export
ExposeAttribute <- function(object, ..., .MODE='rw', .VALUE=FALSE){
	
	# setup exposed arguments
	args <- list(...)
	if( length(args) ){
		# use the same mode for all attributes
		if( isString(.MODE) == 1L )
			.MODE <- rep(.MODE, length(args))
		else if( length(.MODE) != length(args) ){
			stop("Argument .MODE must provide an access mode for each argument in `...`.")
		}
		
		if( is.null(names(args)) ) # add names if necessary
			args <- setNames(args, rep('', length(args)))
		un <- names(args)==''
		if( any(!sapply(args[un], isString)) )
			stop("All unnamed argument must be the name of an attribute, i.e. a character string.")
		
		# set attributes that have values if requested
		if( .VALUE ){
			sapply(names(args)[!un], function(x){
				attr(object, x) <<- args[[x]]
			})
		}else{ # or use the values as access permission
			.MODE[!un] <- args[!un]
		}
		#
		
		# store exposed attributes with names as regular expressions
		eargs <- ifelse(un, args, names(args))
		eargs <- as.list(setNames(.MODE, eargs))
		# add ereg start-end
		names(eargs) <- paste('^', names(eargs), '$', sep='')
	}else{
		eargs <- .MODE
	}
	
	# store access rights
	attr(object, '.ExposeAttribute') <- eargs
	class(object) <- c(class(object), 'ExposeAttribute')
	object
}

.getEAmode <- function(x, name, ..., RAW..=FALSE){
	ea <- attr(x, '.ExposeAttribute')
	if( is.null(ea) ) return()
	if( is.character(ea) && !RAW.. )
		ea <- list(`^.*$`=ea)
	if( missing(name) )	return(ea)

	name <- name[name != '.ExposeAttribute']
	# determine access mode
	m <- lapply(names(ea), function(p){
		m <- grep(p, name, value=TRUE)
		setNames(rep(ea[[p]], length(m)), m)
	})
	unlist(m)
	#
}

#' @importFrom utils .DollarNames
#' @S3method .DollarNames ExposeAttribute 
.DollarNames.ExposeAttribute <- function(x, pattern=""){ 
	
	att <- grep(pattern, names(attributes(x)), value=TRUE)
	if( nchar(pattern) > 1 )
		att <- unique(c(substring(pattern, 2), att))
	# filter out based on the access permissions
	mode <- .getEAmode(x, att)
	if( !length(mode) ) return(character())
	mode <- mode[mode != '']
	# add `<-` suffix to write only attributes
	if( length(wonly <- which(mode=='w')) )
		names(mode)[wonly] <- paste(names(mode)[wonly], '<- ')
	
	names(mode)
}

#' @S3method $ ExposeAttribute
`$.ExposeAttribute` <- function(x, name){
	if( is.null(attr(x, name)) )
		stop("Object `", deparse(substitute(x)),"` has no attribute '", name, "'.")
	mode <- .getEAmode(x, name)
	if( !length(mode) ){
		stop("Could not access attribute via `$`: attribute '", name, "' is not exposed. Use `attr(x, '", name, "').")
	}
	if( !any(grepl('r', mode)) ){
		stop("Could not access exposed attribute '", name, "': permission denied [mode='", mode,"'].")
	}
	attr(x, name)
	
}

#' @S3method $<- ExposeAttribute
`$<-.ExposeAttribute` <- function(x, name, value){
	mode <- .getEAmode(x, name)
	if( !length(mode) ){
		stop("Could not write attribute via `$<-`: attribute '", name, "' is not exposed. Use `attr(x, '", name, "') <- value.")
	}
	if( !any(grepl('w', mode)) ){
		stop("Could not write attribute '", name, "': permission denied [mode='", mode,"'].")
	}
	attr(x, name) <- value
	x
}

#' @S3method print ExposeAttribute
print.ExposeAttribute <- function(x, ...){
	# remove EA stuff 
	attr_mode(x) <- NULL
	# call next print method
	print(x, ...)
}

#' \code{attr_mode} and \code{attr_mode<-} get and sets the access mode of 
#' \code{ExposeAttribute} objects.
#' 
#' @param x an \code{ExposeAttribute} object
#' @param value replacement value for mode.
#' It can be \code{NULL} to remove the ExposeAttribute wrapper, 
#' a single character string to define a permission for all atributes
#' (e.g., \code{'rw'} or \code{'r'}), or a list specifying access permission 
#' for specific attributes or classes of attributes defined by regular expressions.
#' For example, \code{list(a='r', b='w', `blabla.*`='rw')} set attribute \code{'a'} 
#' as read-only, attribute \code{'b'} as write-only, all attributes that start with 
#' \code{'blabla'} in read-write access.
#' 
#' @export
#' @rdname ExposeAttribute
attr_mode <- function(x){
	.getEAmode(x, RAW..=TRUE)
}
#' @export
#' @rdname ExposeAttribute
`attr_mode<-` <- function(x, value){
	if( is.null(value) ){
		attr(x, '.ExposeAttribute') <- NULL
		class(x) <- class(x)[!class(x) %in% "ExposeAttribute"]
	}else if( isString(value) ){
		x <- ExposeAttribute(x, .MODE=value)
	}else if( is.list(value) ){
		args <- c(list(x), names(value), list(.MODE=setNames(value, NULL), .VALUE=FALSE))
		x <- do.call('ExposeAttribute', args)
	}else{
		stop("Invalid value: a character string or a list is expected")
	}
	x
}


#' Checking R User 
#' 
#' Tests if the current R user is amongst a given set of users.
#' 
#' @param user the usernames to check for, as a character vector. 
#' 
#' @export
userIs <- function(user){
	setNames(Sys.info()['user'], NULL) %in% user
}

#' Expanding Lists
#' 
#' \code{expand_list} expands a named list with a given set of default items,
#' if these are not already in the list, partially matching their names.  
#' 
#' @param x input list
#' @param ... extra named arguments defining the default items.
#' A list of default values can also be passed as a a single unnamed argument.
#' @param .exact logical that indicates if the names in \code{x} should be 
#' partially matched against the defaults.
#' @param .names logical that only used when \code{.exact=FALSE} and indicates
#' that the names of items in \code{x} that partially match some defaults should
#' be expanded in the returned list.
#' 
#' @return a list  
#' 
#' @export 
#' @examples 
#' 
#' expand_list(list(a=1, b=2), c=3)
#' expand_list(list(a=1, b=2, c=4), c=3)
#' # with a list
#' expand_list(list(a=1, b=2), list(c=3, d=10))
#' # no partial match
#' expand_list(list(a=1, b=2, c=5), cd=3)
#' # partial match with names expanded
#' expand_list(list(a=1, b=2, c=5), cd=3, .exact=FALSE)
#' # partial match without expanding names
#' expand_list(list(a=1, b=2, c=5), cd=3, .exact=FALSE, .names=FALSE)
#' 
#' # works also inside a function to expand a call with default arguments
#' f <- function(...){
#' 	cl  <- match.call()
#' 	expand_list(cl, list(a=3, b=4), .exact=FALSE)
#' }
#' f()
#' f(c=1)
#' f(a=2)
#' f(c=1, a=2)
#' 
expand_list <- function(x, ..., .exact=TRUE, .names=!.exact){
	
	# extract defaults from ... arguments
	defaults <- list(...)
	if( length(defaults) == 1L && is.null(names(defaults)) ){
		defaults <- defaults[[1L]]
	}
	# early exit if no defaults
	if( !length(defaults) ) return(x)
	
	# match names from x in defaults
	x_ex <- x
	if( !.exact ){
		i <- pmatch(names(x), names(defaults))
		# first expand names if necessary
		if( length(w <- which(!is.na(i))) ){
			names(x_ex)[w] <- names(defaults)[i[w]]
			# apply to as well if necessary
			if( .names ) names(x)[w] <- names(defaults)[i[w]]
		}
	}
	
	# expand list
	i <- match(names(defaults), names(x_ex))
	if( length(w <- which(is.na(i))) ){
		n <- names(defaults)[w]
		lapply(n, function(m){
			if( is.null(defaults[[m]]) ) x[m] <<- list(NULL) 
			else x[[m]] <<- defaults[[m]]
		})
	}
	
	x
}

#' \code{expand_dots} expands the \code{...} arguments of the function
#' in which it is called with default values, using \code{expand_list}.
#' It can \strong{only} be called from inside a function.
#' 
#' @param .exclude optional character vector of argument names to exclude 
#' from expansion. 
#'
#' @export
#' @rdname expand_list
#' 
#' @examples
#' # expanding dot arguments
#' 
#' f <- function(...){ 
#' 	expand_dots(list(a=2, bcd='a', xxx=20), .exclude='xxx') 
#' }
#' 
#' # add default value for all arguments 
#' f()
#' # add default value for `bcd` only
#' f(a=10)
#' # expand names
#' f(a=10, b=4)
#' 
expand_dots <- function(..., .exclude=NULL){
	
	dotsCall <- as.list(eval(quote(substitute(list(...))), sys.parent()))
	if( length(dotsCall) >= 1L ) dotsCall <- dotsCall[-1L]
	
	# extract defaults from ... arguments
	defaults <- list(...)
	if( length(defaults) == 1L && is.null(names(defaults)) ){
		defaults <- defaults[[1L]]
	}
	if( length(defaults) ){
		excl <- names(allFormals(sys.function(sys.parent())))
		if( !is.null(.exclude) ) excl <- c(excl, .exclude)
		defaults <- defaults[!names(defaults) %in% excl]
		dotsCall <- expand_list(dotsCall, defaults, .exact=FALSE)
	}
	#
	
	# return expanded dot args
	dotsCall
}

#' Check Environment Variables
#' 
#' Tells if some environment variable(s) are defined.
#' 
#' @param x environment variable name, as a character vector.
#' 
#' @export
#' @examples 
#' 
#' hasEnvar('_R_CHECK_TIMINGS_')
#' hasEnvar('ABCD')
#' 
hasEnvar <- function(x){
	is.na(Sys.getenv(x, unset = NA, names = FALSE))
}
