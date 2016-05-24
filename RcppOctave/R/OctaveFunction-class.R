# Wrapper class for Octave functions
# 
# Author: Renaud Gaujoux
# Creation: 05 Nov 2011
###############################################################################

#' Wrapping and Defining Octave Functions from R
#' 
#' @slot name name of the wrapped Octave function
#' 
setClass("OctaveFunction", contains="function"
	, representation(name='character')	
)

#' The function \code{OctaveFunction} is a constructor/factory method for 
#' \code{OctaveFunction} objects, which wrap calls to Octave functions into 
#' plain R functions.
#' 
#' \code{OctaveFunction} objects can be created from existing Octave function
#' using their name, or directly from their Octave implementation.
#' In this case, the Octave code is parsed to extract and use the name of the first 
#' function defined therein.
#'  
#' @param fun the name of an existing Octave function or, Octave code that 
#' defines a function.
#' @param check logical that indicates if the existence of the Octave function 
#' should be checked.
#' If function does not exist then, an error or a warning is thrown if \code{check=TRUE}
#' or \code{check=FALSE} respectively.
#' The existence check can be completly disabled with \code{check=NA}.   
#'
#' @rdname OctaveFunction-class
#' @aliases show,OctaveFunction-method
#' 
#' @export
#' @examples
#' 
#' osvd <- OctaveFunction('svd')
#' osvd
#' osvd(matrix(1:9,3))
#' 
#' orand <- OctaveFunction('rand')
#' orand()
#' orand(2)
#' orand(2, 3)
#' 
#' # From source code
#' myfun <- OctaveFunction('function [Y] = somefun(x)
#' 	Y = x * x;
#' 	end
#' ')
#' myfun
#' myfun(10)
#' 
OctaveFunction <- function(fun, check=TRUE){
	
	e <- new.env()
	# use name as text if detected as an octave function defintion
	if( any(d <- is_mdef(fun)) ){
		o_source(text=fun)
		fun <- names(d)[which(d)[1L]]
	}
	
	if( !is_NA(check) && !o_exist(fun) ){
		if( isTRUE(check) ){
			stop("Could not create OctaveFunction object: Octave function '", fun, "' does not exist")
		}else if( isFALSE(check) ){
			warning("Octave function '", fun, "' does not exist: calling the created OctaveFunction object might throw an error.")
		}
	}
	
	# define wrapper function
	f <- evalq({
		.NAME <- fun
		function(...){
			.CallOctave(.NAME, ...)		
		}		
	}, e)
	new('OctaveFunction', f, name=fun)
}


is_mdef <- function(x){
	m <- str_match(x, "((^)|([\n;]))\\s*function\\s*(((\\[([^]]*)\\])|([^=( \t]*))\\s*=)?\\s*([^(]+)")
	setNames(!is.na(m[,1]), m[,10])
}

#' @noRd
#' @export
setMethod('show', 'OctaveFunction', function(object){
	cat("<OctaveFunction::`", object@name, "`>\n", sep='')	
})

#' M Files
#' 
#' \code{as.mfile} converts source code or .m filenames into real paths to .m files
#' that can be sourced with \code{\link[RcppOctave]{o_source}}.
#' 
#' @param ... specification of a .m files as character arguments.
#' The elements of the vector can be either file paths or plain Octave/Matlab code, 
#' which are then written to disk in -- temporary -- .m files. 
#' Note that the paths do not need to correspond to existing files.
#' @inheritParams base::tempfile
#' @param dir existing directory where to write the .m files generated from 
#' the plain code elements of \var{x}.
#' 
#' @export
#' @rdname mfiles
#' @examples 
#' 
#' f <- as.mfile('test.m')
#' f
#' 
#' # detected code elements are written into temporary files
#' f <- as.mfile('test.m', "function [y] = myfun()
#' y = 1;
#' end
#' ")
#' 
#' \dontrun{
#' file.show(f[2])
#' }
#' 
#' # remove all files
#' unlink(f)
#' 
as.mfile <- function(..., pattern='mfile_', dir=tempdir()){
	
	# get args
	x <- unlist(list(...))
	if( !is.character(x) )
		stop("All arguments must be character strings")
	
	# detect type of input
	isfile <- !is_mdef(x)
	# add names if needed
	if( is.null(names(x)) ) names(x) <- rep('', length(x))
	
	code <- x[!isfile]
	if( length(code) ){
        
        in_package <- FALSE
        if( missing(dir) && !is.null(ns <- getLoadingNamespace()) ){
	            in_package <- TRUE
	            dir <- packagePath('m-files', package=ns)
        }
        
		x[!isfile] <- mapply(function(f, x){
					
			# create directory if it does not exist
			if( !file.exists(dir) )	dir.create(dir, recursive=TRUE)
			
			# build file path
			ofile <- f
			if( nchar(f) ) f <- file.path(dir, f)
			else f <- tempfile(pattern, tmpdir=dir)
			f <- str_c(f, '.m')
			# write file
			cat(x, file=f)
			
			# return filepath
			if( in_package ) ofile else f
		}, names(code), code)
	}
	x
}

#' Path to Package M-files Standard Location
#' 
#' \code{system.mfile} returns paths to .m files installed with packages.
#' 
#' \code{system.mfile} is a shortcut for:
#' \samp{system.file('m-files', ..., package = package)}
#' As such it returns empty strings if the requested file does not exist.
#' If no arguments besides \code{package} are passed, it returns the full 
#' path to the package's sub-directory \emph{m-files/} -- if it exists.
#' 
#' @inheritParams base::system.file
#' @param ... arguments passed to \code{\link{system.file}}.
#'  
#' @export
system.mfile <- function(..., package = 'base'){
    system.file('m-files', ..., package = package)
}
