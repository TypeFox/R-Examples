# Functions related to Octave code evaluation
# 
# Author: Renaud Gaujoux
# Creation: 04 Nov 2011
###############################################################################

#' Evaluate an Octave Expression
#' 
#' Evaluates an Octave expression in the current embedded Octave 
#' session. The variables assigned in the expression are available for 
#' subsequent \code{o_eval} calls.
#' 
#' @param ... The Octave expression(s) to evaluate, as a character string.
#' @param CATCH The Octave expression(s) to evaluate if the evaluation(s) of 
#' \code{...} fails. See section \emph{Octave Documentation} for more details.
#' @param unlist a logical that specifies it single variables should be 
#' returned as a single value (default), or as a list.
#'  
#' @return the result of the evaluation
#' 
#' @templateVar name evalin
#' @template OctaveDoc
#' 
#' @export 
#' @examples 
#'
#' # assign some variable 
#' o_eval("a=10")
#' 
#' # retrieve its value in a subsequent call
#' o_eval("a")
#' \dontshow{ stopifnot( identical(o_eval("a"), 10) ) }
#' o_get('a')
#' 
#' # use its value
#' o_eval("b = a^2")
#' \dontshow{ stopifnot( identical(o_eval("b = a^2"), 100) ) }
#' 
#' # multiple expression can be evaluated
#' o_eval(a="10^3", singular="svd(rand(4,4))", random="rand(10, 1)")
#' # or from a list
#' l <- list(a="10^3", singular="svd(rand(4,4))", random="rand(10, 1)")
#' o_eval(l)
#' 
#' # if the evaluation fails then an error is thrown
#' \dontrun{ o_eval("a=svd()") }
#' 
#' # except if argument CATCH is provided 
#' o_eval("a=svd()", CATCH="a=2")
#' \dontshow{ stopifnot( identical(o_eval("a"), 2) ) }
#' 
#' 
o_eval <- function(..., CATCH, unlist=TRUE){
	
	# retrieve the expression to evaluate
	dots <- list(...)
	
	# convert single character vector into a list of evaluations
	if( length(dots) == 1 ){
		if( is.list(dots[[1]]) )
			dots <- dots[[1]]
		else if( is.character(dots[[1]]) ) 
			dots <- as.list(dots[[1]])	
	}
	
	res <- 
			if( missing(CATCH) )
				sapply(dots, .CallOctave, .NAME='eval', simplify=FALSE)
			else 
				sapply(dots, function(x){ .CallOctave('eval', x, CATCH) }, simplify=FALSE)		
	
	# unlist the result if requested
	if( unlist && length(res) == 1 ) res[[1]]
	else res
	
}

#' Assign/Get Octave Variables
#' 
#' \code{o_assign} assigns a variable in Octave.
#' \code{o_assignin} is an alias for \code{o_assign}.
#' 
#' \code{o_assign} assigns the variables using the arguments' names if present.
#' Variables can also be specified as a single named list or environment.
#' Single variable assignments can also be specified as \code{o_assign('a', 10)}.
#' See \emph{Examples} for more details. 

#' 
#' @param ... variables to assign in Octave global context for \code{o_assign}
#' , or object names to retrieve from Octave for \code{o_get}.
#' 
#' @return \code{o_assign} returns invisibly the names of the assigned variables. 
#' 
#' @templateVar name assignin
#' @template OctaveDoc
#' @importFrom stats setNames
#' 
#' @export 
#' @examples 
#' 
#' \dontshow{o_clear()}
#' 
#' ## directly assign variables 
#' o_assign(a=1, b=2, c=matrix(1:9, 3))
#' # retrieve their values
#' o_get()
#' \dontshow{ stopifnot( identical(o_get(), list(a=1, b=2, c=matrix(1:9, 3))) ) }
#' 
#' ## assign a variable for each element in a list
#' x <- list(a=10, b=20, c=matrix(101:109, 3))
#' o_assign(x)
#' o_get()
#' \dontshow{ stopifnot( identical(o_get(), x) ) }
#' 
#' ## assign the content of an environment
#' e <- list2env(setNames(x, paste('env', names(x), sep='_')))
#' o_assign(e)
#' o_get(pattern="^env_")
#' 
#' 
o_assign <- function(...){
	
	# check that all arguments are named
	VALUES <- list(...)
	
	# convert unnamed list or environment into a "list" 
	if( is.null(names(VALUES)) ){
		
		if( length(VALUES) == 1 ){
						
			if( is.list(VALUES[[1]]) || is.environment(VALUES[[1]]) )
				VALUES <- VALUES[[1]]
			else 
				VALUES <- as.list(VALUES[[1]])			
			
		}else if( length(VALUES) == 2 ){
			n <- eval(VALUES[[1]])
			VALUES <- setNames(VALUES[2], as.character(substitute(n)))
		}else
			stop("Invalid assignment command: expecting a single named list or environment, 2 unnamed arguments, or named arguments.")
	}	
	
	VARNAMES <- if( is.environment(VALUES) ) ls(VALUES) else names(VALUES)
	if( is.null(VARNAMES) || any(VARNAMES=='') )
		stop("All arguments must be named to assign.")
	
	CONTEXT='base'
	sapply(VARNAMES, function(name){
				.CallOctave('assignin', CONTEXT, name, VALUES[[name]])
			}, simplify=FALSE)
	invisible(VARNAMES)
}
#' @export
#' @rdname o_assign 
o_assignin <- o_assign

#' \code{o_get} retrieves objects from Octave.
#' 
#' \code{o_get} fetches Octave variables/functions and possibly rename them on 
#' the fly with the provided argument names when present.
#' Functions are returned as objects of class \code{\linkS4class{OctaveFunction}}, 
#' that can be called subsequently (see the examples).
#' 
#' @note The function \code{o_get} is the equivalent of R \code{\link{get}} 
#' function and is not related to the Octave function \code{get} which returns 
#' graphics' properties.
#' 
#' @param unlist a logical that specifies it single variables should be 
#' returned as a single value (default), or as a list.
#' @param rm.ans a logical that indicates if the automatic Octave variable \code{ans}
#' should be included in the result. Default is not to include it unless otherwise 
#' explicitly specified with this argument, or if it is part of the requested 
#' variables in \code{...}. When present, argument \code{rm.ans} is always honoured.
#' @param pattern regular expression used to filter the requested variable names.
#' Only names matching the pattern are returned. 
#' 
#' @return \code{o_get} returns a list of the retrieved variable/object. If 
#' \code{unlist=TRUE} and a single -- not re-named -- variable/object is requested then
#' only its value is returned.
#' 
#' @rdname o_assign
#' @export
#' @examples 
#' 
#' \dontshow{o_clear();}
#' 
#' # get all currently defined variables
#' o_get()
#' 
#' # by default, the automatic variable `ans` is not returned but might be there
#' # from unstored previous computation 
#' o_eval('svd(rand(3,3))')
#' o_get()
#' o_get(rm.ans=FALSE)
#' 
#' # load some variables
#' x <- list(b=1, c=3, d=matrix(1:9, 3))
#' o_assign(x)
#' 
#' # re-fetch all variables
#' o_get()
#' \dontshow{ stopifnot( identical(o_get(), x) ) }
#' 
#' # only fetch specific variables
#' o_get('b')
#' o_get('b', 'c')
#' # one can rename variables on the fly
#' o_get(a='b', 'c')
#' o_get(c(othername='b'))
#' o_get(c(othername='b', 'c'))
#' 
#' # fetching using a regular expression 
#' o_assign( setNames(1:3, paste("test", 1:3, sep='_')))
#' o_get()
#' o_get(pattern="^test")
#' 
#' # works with functions
#' f <- o_get('svd')
#' f
#' f(matrix(1:9,3))
#' f(matrix(1:9,3), argout=3)
#' 
#' # an error is thrown in the case of multiple matches (the alternatives are shown)
#' o_clear()
#' o_assign(aaa=1, ab=2)
#' try(o_get('a'))
#' 
o_get <- function(..., unlist=TRUE, rm.ans = TRUE
#				, exact=FALSE
				, pattern){
			
	dots <- unlist(list(...))
	vnames <- dots
	#print(dots)
	if( is.null(names(dots)) || length(names(dots)) == 0 ){
		vnames <- 
		if( length(dots) == 0 ){ # no argument: get all variables
			
			# enforce no ans if not otherwise requested
			if( missing(rm.ans) )
				rm.ans <- TRUE
			vnames <- o_who()
			setNames(vnames, vnames)
			
		}
		else if( is.character(dots) )
			dots
		else if( length(dots) == 1 ){			
			# single unnamed argument: use it as a vector of names			
			dots[[1]]
		}else
			stop("Invalid get command with unnamed arguments: expecting no arguments or character vectors.")		
	}
		
	#print(vnames)
	# keep 'ans' if it is part of the request and not otherwise explicitly requested 
	if( missing(rm.ans) && 'ans' %in% vnames )
		rm.ans <- FALSE
	# add names if necessary 
	if( is.null(names(vnames)) ){
		vnames <- setNames(vnames, vnames)
	}else if( missing(unlist) ){
		# disable unlisting if names are specified and not otherwise explicitly requested
		unlist <- FALSE
		noname <- names(vnames) == ''
		names(vnames)[noname] <- vnames[noname]
	}
	# remove the automatic variable 'ans' if requested	
	if( rm.ans )
		vnames <- vnames[vnames != 'ans']
	# subset using the pattern
	if( !missing(pattern) )
		vnames <- vnames[grep(pattern, vnames)]
	
	if( length(vnames) == 0 ) return(list())
	
	# get the value for each name
	res <- lapply(seq_along(vnames), function(i){
		
		name <- vnames[[i]]
		# check for an exact match
		ecode <- o_exist(name)
		if( !ecode ){
			
			onames <- o_completion_matches(name)
			sugg <- 
			if( length(onames) > 0L )
				str_c("\n       Match(es): ", str_wrap(paste(onames, collapse=" "), exdent=18), "\n")
			
			stop("RcppOctave::o_get - Could not find an Octave object named '", name,"'.", sugg, call.=FALSE)
			
#			if( exact )
#				stop("RcppOctave::o_get - Could not find an Octave object named '", name,"'.", call.=FALSE)
#			
#			# check for multiple matches
#			onames <- o_completion_matches(name)
#			if( length(onames) == 0 )
#				stop("RcppOctave::o_get - Could not find an Octave object starting with '", name,"'.", call.=FALSE)
#			
#			if( length(onames) > 1 )
#				stop("RcppOctave::o_get - Multiple Octave objects [", length(onames), "] start with '", name,"'.\n"
#						, "       Matches are: ", stringr::str_wrap(paste(onames, collapse=" "), exdent=20), "\n"
#						, call.=FALSE)
#			name <- onames
		}
		
		# check if `name` is a variable
		if( ecode == 1 ){ # return the value of the variable
			o_eval(name)
		}else # wrap Octave functions into R wrappers
			OctaveFunction(name)
	})
	res <- setNames(res, names(vnames))

	# unlist the result if requested
	if( unlist && length(res) == 1 ) res[[1]]
	else res
}


#' Loading Variables into Octave
#' 
#' Loads variables from a file, a list or an environment.
#' 
#' @param from a list or an environment from where the objects should be loaded
#' @param ... names of the variables to load
#' @param options argument passed to the Octave function \code{load}. See section 
#' \emph{Octave Documentation}.
#'  
#' @templateVar name load
#' @template OctaveDoc
#' 
#' @export
#' @examples
#' 
#' # Loading from a MATLAB/Octave file
#' #o_load
#' 
#' # Loading from an R list
#' o_clear()
#' l <- list(a=1, b=20, c=runif(10), d="this is a string", e=matrix(1:15, 3, 5))
#' o_load(l)
#' 
#' # Loading from an R environment
#' o_load( list2env(l) )
#' 
#' # Partial loading
#' o_clear()
#' o_load(l, a, b, c)
#' o_clear()
#' o_load(list2env(l), d, e)
#' 
o_load <- function(from, ..., options){
	
	# use o_assign for lists and environments
	if( is.list(from) || is.environment(from) ){
		names <- as.character(substitute(list(...)))
		# subset values if requested
		if( length(names) > 1 ){
			names <- names[-1]
			from <- 
			if( is.list(from) ) from[names]
			else mget(names, from)
			
		}
		return(o_assign(from))
	}
	
	# check file existence
	if( !file.exists(from) )
		stop("Could not load variables: file '", from , "' does not exist.")
	
	if( missing(options) ) .CallOctave('load', from, ...)
	else  .CallOctave('load', options, from, ...)
}

#' Deleting Octave Variables
#' 
#' Deletes variables from Octave global context.
#' 
#' @param ... names or pattern of the variables to delete, as character strings.
#' @param all a logical indicating whether all user-defined objects should be 
#' deleted. See section \emph{Octave Documentation} for details.
#' @param options options passed to Octave function \code{clear}. 
#' See section \emph{Octave Documentation}.
#' 
#' @return None
#' 
#' @templateVar name clear
#' @template OctaveDoc
#' 
#' @export
#' @examples
#' \dontshow{ o_clear() }
#' # Assign a variable in Octave
#' o_assign('a', 10)
#' o_who()
#' \dontshow{ identical(o_who(), 'a') }
#' # Clear
#' o_clear()
#' o_who()
#' \dontshow{ identical(o_who(), character()) }
#' 
#' # Assign other variables in Octave
#' .O$a <- 10
#' .O$b <- 100
#' .O$ba <- 1000
#' o_who()
#' o_get()
#' \dontshow{ identical(o_who(), c('a', 'b', 'ba')) }
#' 
#' # Clear variable starting with 'b'
#' o_clear('b*')
#' o_who()
#' \dontshow{ identical(o_who(), 'a') }
#' 
#' 
o_clear <- function(..., all=FALSE, options){
	
	# add option -all
	if( all ){
		if( missing(options) ) options <- '-all'
		else options <- paste(options, '-all')
	}

	if( missing(options) ) .CallOctave('clear', ...)
	else  .CallOctave('clear', options, ...)
	invisible()
}

#' The function \code{o_rm} is an alias to \code{o_clear}.
#' @rdname o_clear
#' @export
o_rm <- o_clear

#' Listing Octave Variables
#' 
#' Lists currently defined variables in Octave global context.
#' 
#' @param ... filtering patterns or extra arguments passed to \code{o_who} 
#' and \code{o_whos}. Only names matching any of the patterns are returned.
#' @param rm.ans a logical that indicates if the automatic Octave variable \code{ans}
#' should be included in the result (\code{FALSE}) or removed (\code{TRUE}).
#' @param options options passed to Octave function \code{who}. 
#' See section \emph{Octave Documentation}.
#' @param unique a logical that indicates whether unique names should be returned.
#' This argument is relevant in the case multiple patterns are specified in 
#' \code{...}.
#' 
#' @return None
#' 
#' @templateVar name who
#' @template OctaveDoc
#' 
#' @family listoct 
#' @export
#' @examples 
#' 
#' \dontshow{ o_clear() }
#' o_who()
#' l <- as.list(setNames(1:10, letters[1:10]))
#' o_assign(l)
#' o_who()
#' \dontshow{  stopifnot( identical(o_who(), names(l)) ) }
#' 
#' prefnames <- paste('pref', letters[1:10], sep='')
#' o_assign( setNames(l, prefnames) )
#' o_who()
#' o_who('pref*')
#' \dontshow{  stopifnot( identical(o_who('pref*'), prefnames) ) }
#' 
#' 
o_who <- function(..., options, rm.ans=FALSE, unique=TRUE){
	
#	# add option -long
#	if( long ){
#		if( missing(options) ) options <- '-long'
#		else options <- paste(options, '-long')
#	}
#	
#	# add option -all
#	if( all ){
#		if( missing(options) ) options <- '-all'
#		else options <- paste(options, '-all')
#	}
	
	res <- 
	if( missing(options) ) .CallOctave('who', ...)
	else  .CallOctave('who', options, ...)

	# remove automatic variable 'ans'
	if( rm.ans ) res <- res[res != 'ans']

	if( unique ) unique(res)
	else res
}

#' Detailed Listing of Octave Variables
#' 
#' The function \code{o_whos} returns a detailed description of the variables
#' defined in the current Octave session.
#' 
#' @inheritParams o_who
#'  
#' @templateVar name whos
#' @template OctaveDoc
#' 
#' @family listoct
#' @export
#' 
#' @examples
#' \dontshow{ o_clear() }
#' 
#' .O$a <- 1
#' .O$b <- 10
#' o_whos()
#' 
#' o_eval("sqrt(2)")
#' o_whos()
#' 
o_whos <- function(..., options, rm.ans=FALSE){

	ofun <- 'whos'
	w <- 
	if( missing(options) ) .CallOctave(ofun, ...)
	else  .CallOctave(ofun, options, ...)
	
	# special case to ensure to have a lists in 'size' and 'nesting' 
	if( length(w$name) == 1L ){
		w$size <- list(w$size)
		w$nesting <- list(w$nesting)
	}
	
	# remove automatic variable 'ans' if necessary
	if( rm.ans ){
		wans <- which(w$name == 'ans')
		if( length(wans) > 0L )
			w <- lapply(w, '[', -wans)
	}
	
	# add S3 class for better formating 
	class(w) <- c("Octave_whos", class(w))
	w
}

#' @S3method print Octave_whos
#' @noRd 
print.Octave_whos <- function(x, ...){
	# reformat some of the data
	x$size <- sapply(x$size, paste, collapse="x")
	x$nesting <- sapply(x$nesting, paste, collapse="")
	x <- as.data.frame(x)
	cat(" <Octave session: ", nrow(x)," object(s)>\n", sep='')
	if( nrow(x) > 0L ) print(x, ..., row.names=FALSE)
	invisible(x)
}

#' Listing Objects from the Current Octave Session
#'  
#' The function \code{o_ls} is an enhanced listing function, which 
#' also lists user-defined functions, as opposed to \code{\link{o_who}} or 
#' \code{\link{o_whos}}, which only list variables.
#' Note that this function works properly on Octave >= 3.6.1, but is known 
#' not to list user-defined functions on Octave 3.4.1 (for some unknown reason the 
#' Octave function \code{completion_matches} does not return the names of 
#' user-defined functions).
#' 
#' @param rm.ans a logical that indicates if the automatic Octave variable \code{ans}
#' should be included in the result. Default (\code{TRUE}) is not to include it.
#' 
#' @param long logical that indicates the desired type of output: if \code{FALSE} 
#' (default) then only the names of the variables are returned (like dispatched 
#' \code{o_who}), otherwise a list with more detailed information about 
#' each variable is returned (like \code{\link{o_whos}}.
#' 
#' @return a character vector or a list depending on the value of argument 
#' \code{long}.
#' 
#' @import digest
#' @family listoct
#' @export
#' @examples 
#' 
#' \dontshow{ o_clear(all=TRUE) }
#' 
#' # only variables
#' o_assign(list(a=1, b=2, c=5))
#' o_ls()
#' # compare with the output of standard Octave functions
#' o_who() # should be the same output
#' o_whos()
#' 
#' # variables and user-defined functions
#' o_clear(all=TRUE) # first clear Octave session
#' o_source(system.file('scripts/ex_source.m', package='RcppOctave'))
#' o_ls()
#' o_ls(long=TRUE)
#' # compare with the output of standard Octave functions
#' o_who()
#' o_whos()
#' 
o_ls <- local({
			
	.cache <- NULL
	.cacheMD5 <- ""
	function(long=FALSE, rm.ans=TRUE){
	
		# get all variables
		if( long ){
			var <- o_whos(rm.ans=rm.ans)
			varnames <- var$name
		}else{
			var <- o_who(rm.ans=rm.ans)
			varnames <- var
		}
		# get all objects
		ol <- o_completion_matches("")
		# remove already known variables
		ol <- ol[ !is.element(ol, varnames) ]
		if( {hash <- digest(ol)} != .cacheMD5 ){
			.cacheMD5 <<- hash
			# filter for functions
			ecode <- sapply(ol, o_exist)
			.cache <<- ol[ecode == 103]
		}
		
		if( long ){
			l <- list()
			sapply(.cache, function(x){
				var$name <<- c(var$name, x)				
				var$size <<- c(var$size, list(as.character(NA)))
				var$bytes <<- c(var$bytes, as.numeric(NA))
				var$class <<- c(var$class, "function")
				var$global <<- c(var$global, TRUE)
				var$sparse <<- c(var$sparse, NA)
				var$complex <<- c(var$complex, NA)
				var$nesting <<- c(var$nesting, list(list(`function`="", level=1)))
				var$persistent <<- c(var$persistent, NA)
			})
			var
		}else c(var, .cache)
	}
})	
