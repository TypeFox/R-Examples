# Class Octave: direct interface to Octave functions and base context
# 
# Author: Renaud Gaujoux
# Creation: 05 Nov 2011
###############################################################################


#' Class Octave: Seamless Access to Octave Functions and Variables
#' 
#' This class provides a direct interface to Octave functions and base context
#' (i.e. the scope where user objects are defined).
#'  
#' 
#' @import methods
#' @keywords internal
setClass("Octave", contains='character')

#' Direct Interface to Octave 
#' 
#' \code{RcppOctave} provides a simple interface to Octave via the  
#' object \code{.O}, an instance of class \code{Octave}, that allows for direct access
#' to Octave functions and variables using calls such as: \code{.O$svd(matrix(1:9,3))}.
#' 
#' 
#' @rdname OctaveInterface
#' @aliases show,Octave-method
#' @format \code{.O} is an object of class \code{\linkS4class{Octave}}.
#' @export
#' @examples
#' .O
#' # assign/get Octave variables
#' .O$a <- 10
#' .O$a
#' 
#' # call Octave functions
#' .O$help()
#' .O$svd(matrix(runif(9), 3))
#' 
.O <- new("Octave")

#' @noRd 
#' @export
setMethod('show', 'Octave',
		function(object){
			cat(" <Octave Interface>\n")
			cat(" - Use `$x` to call Octave function or get variable x.\n")
			cat(" - Use `$x <- val` to assign a value val to the Octave variable x.\n")
		}
)

setGeneric('.DollarNames', package='utils')

o_completion_matches <- function(pattern = ""){
	# remove leading "^"	
	opattern <- sub("^[\\^]?(.*)", "\\1", pattern)
	grep(pattern, .CallOctave('completion_matches', opattern), value=TRUE)	
}

#' Checking Octave Variables
#' 
#' Checks if an Octave object of a given name exists,
#' using the Octave function \code{exist}.
#'
#' @templateVar name exist
#' @template OctaveDoc
#'
#' @param NAME name to check existence.
#' @param ... extra parameters passed to the Octave function 
#' \code{exist}.
#'
#' @export  
o_exist <- function(NAME, ...){

	ecode <- .CallOctave('exist', NAME, ...)
	if( !length(ecode) ){
		stop("Unexpected empty result from Octave function 'exist' [Octave version: ", o_version(), "].")
	}
		
	if( o_version("3.2.4") == 0 ){ # special handling in 3.2.4
		if( ecode == 0 && NAME %in% o_completion_matches(NAME) ){
			ecode <- 1
		}
	}
	ecode
	
}

#' Auto-completion for Octave Interface Object
#' 
#' Auto-completion for \code{\linkS4class{Octave}} objects
#' 
#' @inheritParams utils::.DollarNames
#' @S3method .DollarNames Octave
#' @rdname autocomplete
#' @keywords internal
.DollarNames.Octave <- function(x, pattern = "") o_completion_matches(pattern)

#' @export
#' @rdname autocomplete
setMethod('.DollarNames', 'Octave', .DollarNames.Octave)

#' The method \code{$} provides a direct way of calling Octave functions or 
#' retrieving variables from Octave base context, via e.g. \code{.O$svd(x)} 
#' or \code{.O$a}.
#' It is equivalent to \code{o_get(name)}. 
#'  
#' 
#' @param x an \code{OctaveInterface} object. Essentially used with \code{x = .O}.
#' @param name name of the Octave object to retrieve or assign.
#' 
#' @rdname OctaveInterface
#' @seealso \code{\link{o_get}}
#' @export
setMethod('$', 'Octave', function(x, name)	o_get(name))


#' The method \code{$<-} allow to directly assign/set Octave variables via e.g.
#' \code{.O$a <- 10}. 
#' 
#' @param value value to assign to the Octave object.
#' @rdname OctaveInterface
#' @export 
setReplaceMethod('$', 'Octave',
	function(x, name, value){
		# remove variable if value is directly NULL
		if( is.null(value) ){
			o_clear(name)			
		}else{
			# force evaluation now
			value <- force(value)
			# assign result of evaluation
			o_assign(name, value)
		}
		x
	}
)

#' The method \code{[[} provides an alternative way of retrieving Octave objects,
#' and is equivalent to \code{o_get(name)}.
#' 
#' @param i name of the Octave object to retrieve.
#' @param exact logical not used.
#'  
#' @rdname OctaveInterface
#' @seealso \code{\link{o_get}}
#' @export
setMethod('[[', 'Octave', function(x, i, exact=TRUE)	o_get(i))


