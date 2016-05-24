# R functions to interact with an embedded Octave instance
# 
# Author: "Renaud Gaujoux"
# Creation: 26 Oct 2011
###############################################################################

#' @include utils.R
NULL

#' Calling an Octave Function
#' 
#' \code{.CallOctave} calls an Octave function and returns its value.
#' 
#' @param .NAME an Octave function name. The function must be a valid function 
#' name in the current Octave session.
#' @param ... arguments passed to the Octave function
#' @param argout the number of output values, or a vector of names to use as output
#' variable names. The names are directly used and applied to the result list in 
#' their original order.
#'  
#' The default value \code{argout=-1} returns:
#' \itemize{
#' \item all output values when their number can be determined. This would typically  
#' be the case for functions defined in .m files. Please do read section 
#' \emph{Details} for considerations about the functions that use varargout. 
#' \item only the first one otherwise.
#' }
#' @param unlist a logical that specifies if an output list of length one 
#' should be simplified and returned as a single value or kept as a list.
#' The default is to unlist unless output names were passed in \code{argout}.
#' @param buffer.std logical that indicates if Octave stdout and/or stderr should be buffered.
#' If \code{TRUE} output/errors/warnings are all displayed at the end of the computation.
#' If \code{FALSE} they are directly displayed as they come.
#' It is also possible to selectively buffer either one of stdout or stderr, via 
#' the following integer codes:
#' \itemize{
#' \item \code{0}: no buffering; 
#' \item \code{1} or \code{-2}: only stdout is buffered;
#' \item \code{2} or \code{-1}: only stderr is buffered;
#' \item \code{4} or \code{-3}: only warnings are buffered;
#' \item \code{7}: all messages are buffered.
#' }
#' 
#' Note that warnings are handle slightly differently than other messages, 
#' as they are never output directly, except when \code{buffer.std = 0}.
#' 
#' @param verbose logical that toggles verbosity (i.e. debug) messages.
#' If \code{NULL}, then the current verbosity level is used 
#' (see \code{\link{octave_verbose}}).
#' 
#' @return the value returned by the Octave function -- converted into standard 
#' R objects.
#' 
#' @export
#' @examples 
#' 
#' # data matrix
#' x <- matrix(1:9, 3)
#' 
#' # call Octave function 'svd': equivalent to [S] = svd(x). See o_help(svd)
#' .CallOctave('svd', x)
#' 
#' # call Octave function 'svd' asking for 3 output values: [U, S, V] = svd(x)  
#' .CallOctave('svd', x, argout=3)
#' 
#' # call Octave function 'svd' asking for 3 named output values: [U, S, V] = svd(x)
#' .CallOctave('svd', x, argout=c('U', 'S', 'V'))
#' 
.CallOctave <- function(.NAME, ..., argout=-1, unlist=!is.character(argout), buffer.std = -1L, verbose = NULL){
	
    if( !is.null(verbose) ){
        ov <- octave_verbose(verbose)
        on.exit( octave_verbose(ov), add = TRUE)
    }
    if( isTRUE(buffer.std) ) buffer.std <- 7L
    res <- .Call("octave_feval", .NAME, list(...), argout, unlist, buffer.std, PACKAGE='RcppOctave')
	if( identical(argout, 0) || identical(argout, 0L) )	invisible()
	else if( is.null(res) && argout <= 1L ) invisible()
	else res
}

#' Low-level Function Interfacing with Octave
#' 
#' \code{octave_start} Initialize an Octave session.
#' 
#' @param verbose logical that toggle verbosity.
#' In \code{octave_start}, it is the value used as the initial global verbosity state. 
#' If \code{TRUE} all calls and conversions between R and Octave produce diagnostic messages.
#' @param warnings logical that indicates if Octave startup warnings
#' @param force logical that indicates if Octave session should be reinitialised, 
#' even if one was previously started (not meant to be used by end-users).  
#' should be shown.
#' 
#' @rdname octave-ll
#' @export
octave_start <- local({
    .Initialised <- FALSE
    function(verbose=FALSE, warnings = FALSE, force = FALSE){
        res <- FALSE
        if( !.Initialised || force ){
	        res <- .Call("octave_start", verbose, warnings, PACKAGE='RcppOctave')
            .Initialised <<- TRUE
        }
        res
    }
})
#' \code{octave_end} clears and terminates the current Octave session.
#' 
#' @rdname octave-ll
#' @export
octave_end <- function(verbose = getOption('verbose')){
	.Call("octave_end", verbose, PACKAGE='RcppOctave')
}

#' \code{octave_verbose} toggles the verbosity of RcppOctave calls: messages tracks 
#' any function call, or conversion of objects between R and Octave 
#' (e.g. arguments and results).
#' 
#' @param value logical value to toggle verbosity
#' 
#' @rdname octave-ll
#' @export
octave_verbose <- function(value){
    value <- if( !missing(value) ) value
	res <- .Call("octave_verbose", value, PACKAGE='RcppOctave')
    if( !is.null(value) ) invisible(res) else res
}

#' Octave Utils: octave-config
#' 
#' Retrieves Octave configuration variables using \code{octave-config}. 
#' 
#' \code{octave_config} uses the \code{octave-config} utility binary shipped with 
#' Octave to query details about the local Octave installation.
#' Failure to retrieve such information is generally due to the binary
#' not being found.
#' By default, it is looked up in the \code{bin/} sub-directory of the path 
#' returned by \code{\link{Octave.home}()}.
#' 
#' @param varname Name (as a character string) of the Octave configuration 
#' variable to retrieve. It is used in following system call 
#' \samp{octave-config -p <varname>}.
#' This function is vectorised and returns a character vector of the same length
#' as its argument.
#' @param verbose logical that toggles verbose messages.
#' @param warn logical that indicates if a warning should be thrown when a 
#' variable is returned empty, which generally means that \code{x} is not a valid 
#' config variable name.
#' @param mustWork logical that indicates if an error should be thrown if failing 
#' to load Octave configuration.
#' @param exec name of the executable to query
#' @param bindir path to Octave bin/ sub-directory where to look for \code{octave-config}.
#' If \code{NULL} or \code{NA}, then the system PATH is looked up.
#'  
#' @family Octave.info
#' @export
#' @examples
#' octave_config('VERSION') 
#' octave_config('BINDIR')
#' 
octave_config <- function(varname, verbose=FALSE, warn=TRUE, mustWork = TRUE, exec = c('octave-config', 'mkoctfile'), bindir = Octave.home('bin')){

    # use custom BINDIR if requested
    octave_config_cmd <- match.arg(exec)
    if( !is.null(bindir) && !is_NA(bindir) ){
        if( verbose ) message("# Using Octave BINDIR '", bindir, "'")
        octave_config_cmd <- file.path(normalizePath(bindir), octave_config_cmd) 
    }

	tryCatch({
        sapply(varname, function(x){
                
            # run octave-config command
            if( verbose ) message("# Loading Octave config variable '", x, "' ... ", appendLF=FALSE)
            cmd <- paste('"', octave_config_cmd, '" -p ', x, sep = '')
    		res <- system_call(cmd)
                       
            # check result
    		if( res == '' ){
    			if( verbose ) message("WARNING")
    			if( warn ) warning("Octave config variable '", x, "' is empty")
    			return(res)
    		}
    		if( verbose ) message('OK')
    		res
    	})
    }
    , error = function(e){
        if( mustWork ) stop(e)
        if( warn ) warning("Failed loading Octave configuration: ", e)
        NULL
        }
    )
}

#' \code{octave_modules} add the Octave modules that ship with RcppOctave to 
#' Octave loading path.
#' 
#' @rdname octave-ll
#' @export
octave_modules <- function(verbose=getOption('verbose')){
	
	path <- Octave.info('modules')
    if( verbose )
		message("Loading Octave modules for ", packageName()
				, " from '", path, "'");
	o_addpath(path)
}

#' Loading Example M-files
#' 
#' Source an example M-file in the sub-directory \dQuote{scripts/} of RcppOctave
#' installation. 
#' 
#' @param file filename of the example script to source. If missing, the function 
#' lists all the m-files from the \dQuote{scripts/} sub-directory. 
#' 
#' @export
#' @examples 
#' 
#' sourceExamples()
#' sourceExamples('ex_source.m')
#' 
sourceExamples <- function(file){
	if( missing(file) ){
		list.files(packagePath('scripts'), pattern="\\.m$")
	}else{# source script
		o_source(packagePath('scripts', file))
	}
	
}


