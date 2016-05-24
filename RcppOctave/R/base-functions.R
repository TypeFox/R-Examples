# Base wrapper functions to some Octave functions
# 
# Author: "Renaud Gaujoux"
# Creation: 27 Oct 2011
###############################################################################

#' Sourcing Octave/Matlab Files
#' 
#' This function sources an Octave file within the current Octave session.
#' The loaded functions are accessible by subsequent calls of 
#' \code{\link{.CallOctave}}.
#' 
#' @param file the path to the Octave/Matlab source file -- typically with 
#' extension ".m".
#' @param text a character vector containing \emph{Octave} statements, that are 
#' concatenated in a temporary file, which is then sourced.
#' This argument typically enables the evaluation of multiple statements, as 
#' opposed to single statement evaluation performed by \code{\link{o_eval}}.
#' @param sep single character string added as suffix to each element of 
#' \code{text}. The concatenation of all suffixed element should form a valid 
#' \emph{Octave} block.
#' 
#' @templateVar name source
#' @template OctaveDoc
#' 
#' @return None
#' @family Octave_files 
#'  
#' @export
#' 
#' @examples 
#' 
#' \dontshow{ o_clear() }
#' 
#' # source file
#' mfile <- system.file("scripts/ex_source.m", package='RcppOctave')
#' o_source(mfile)
#' 
#' # pass multiple statements
#' o_source(text="a=1;b=3;c=randn(1,5);")
#' o_get('a','b','c')
#' 
#' # also works with a character vector of statements
#' o_source(text=c("a=10;b=30;", "c=randn(1,5)", "d=4"))
#' o_get('a','b','c', 'd')
#' 
o_source <- function(file = "", text=NULL, sep=";\n"){
	
	# create temporary file with provided text
	if( !missing(text) ){
		text <- paste(text, sep, sep='')
		mfile <- tempfile("mfile_")
		on.exit(file.remove(mfile))
		cat(text, file=mfile, sep='')
		file <- mfile
	}
	
	if( !file.exists(file) )
		stop("File `", file, "` does not exist.")
	
	.CallOctave('source', file)
	invisible()
}

#' Manipulating Octave Search Path 
#' 
#' Adds a directory at the beginning of Octave search path.
#' 
#' The .oct files present in directories from the search path are looked up 
#' when an object or function is requested but not loaded in the current session.
#' The files are watched and automatically reloaded in case modification.
#' 
#' @param DIR1 path specification to add to Octave search path. 
#' See section \emph{Octave Documentation}.
#' @param ... other path specifications
#' @param OPTION option that specifies how the path should be added. 
#' Possible values are: \code{'-begin', 0, '-end', 1}.
#' See section \emph{Octave Documentation}.
#' 
#' @return returns invisibly the old value of search path.
#' @family Octave_files
#' @export 
#' 
#' @templateVar name addpath
#' @template OctaveDoc
#' 
#' @examples
#' 
#' # call an undefined function
#' try(.CallOctave('fun1'))
#' 
#' # add to the path a directory with a .oct file that contains a definition for 'fun1'
#' o_addpath(system.file('scripts', package='RcppOctave'))
#' 
#' # re-call the function
#' #.CallOctave('fun1')
#' 
#' # change the .oct file
#' 
#' 
o_addpath <- function(DIR1, ..., OPTION='-begin'){
	invisible(.CallOctave('addpath', DIR1, ..., OPTION))
}

o_rmpath <- function(DIR1, ...){
	invisible(.CallOctave('rmpath', DIR1, ...))
}

#' Get Octave Version
#' 
#' Returns the version of Octave currently used by \code{RcppOctave}.
#' 
#' @param version optional reference version to compare with.
#' @return Octave version as a single character string or the result of 
#' \code{\link[utils]{compareVersion}} if argument \code{version} is provided.
#' 
#' @importFrom utils compareVersion
#' @export
#' @examples
#' 
#' o_version()
#' o_version("3.6.2")
#' o_version("3.4")
#' 
o_version <- function(version){
	v <- .CallOctave('version')
	if( !missing(version) ) compareVersion(v, version)
	else v
}

#' Embedded Octave Configuration and Installation Information
#' 
#' Retrieves configuration and installation information about 
#' the embedded instance of Octave currently used by \pkg{RcppOctave}. 
#' 
#' @param var name of the configuration variable to retrieve.
#' 
#' @note 
#' In Octave 4.0.0., variables \code{CC_VERSION} and \code{CXX_VERSION} 
#' were renamed \code{GCC_VERSION} and \code{GXX_VERSION} respectively.
#' \code{o_config_info} supports both values, querying the right ones
#' depending on the linked Octave version.  
#'  
#' See \url{http://www.gnu.org/software/octave/NEWS-4.0.html}
#' 
#' @family Octave.info
#' @export 
#' @examples
#' 
#' o_config_info()
#' o_config_info('USE_64_BIT_IDX_T')
#' 
o_config_info <- function(var = c('CC', 'CC_VERSION', 'FC')){
    
    # fix for Octave version >= 4.0.0: CC_VERSION is now GCC_VERSION
    if( o_version("4.0.0") >= 0 ){
        var0 <- var
        var[ var == "CC_VERSION" ] <- "GCC_VERSION"
        var[ var == "CXX_VERSION" ] <- "GXX_VERSION"
        var[ var == "FC" ] <- "F77"
        names(var) <- var0
    }
    sapply(var, .CallOctave, .NAME = 'octave_config_info')
}

o_builtin <- function(name, ...){
    .CallOctave('builtin', name, ...)
} 

#' Accessing Octave Help and Documentation Pages
#' 
#' \code{o_help} retrieves the Octave help page associated with a given symbol.
#' By default the page is printed out, but may also be silently retrieved or 
#' formatted for direct inclusion in R documentation files (i.e. Rd files). 
#' 
#' @param NAME Octave symbol (e.g. command, function, operator) passed to Octave 
#' function \code{help} to retrieve the related documentation.
#' @param character.only a logical indicating whether \code{NAME} can be
#' assumed to be a character string (\code{TRUE}) or should be substituted
#' with \code{\link{substitute}} before using them (default).
#' @param show logical that specifies if the help page should be shown using the 
#' as R documentation file (default), e.g. using a pager, or only returned as a 
#' single string. Note that when \code{show=TRUE}, the string is still returned 
#' but invisibly.
#' @param format a specification of the output format.
#' If \code{TRUE} or \code{'rd'}, the result is Rd code that wraps the Octave documentation string 
#' and is suitable for inclusion into Rd files.
#' If one of the strings \code{'txt', 'latex'} or \code{'HTML'}, then the result is 
#' formated using the corresponding Rd conversion function from the \pkg{tools} package 
#' \code{\link{Rd2txt}}, \code{\link{Rd2latex}} or \code{\link{Rd2HTML}}.
#' 
#' @return this function is usually called for its side effect of printing the 
#' help page on standard output (argument \code{show=TRUE}), but it invisibly
#' returns the help page as a single character string.
#' 
#' @templateVar name help
#' @template OctaveDoc
#' 
#' @import tools
#' @export
#' @examples 
#' 
#' \dontshow{
#' 	if( interactive() ){
#' 		o_help <- function(..., show=FALSE){
#' 			RcppOctave::o_help(..., show=show)
#' 		}
#' 	}
#' }
#' 
#' o_help(print)
#' o_help(rand)
#' # or equivalently 
#' o_help('rand')
#' 
#' # to include in Rd files, use argument rd=TRUE in an \Sexpr:
#' \dontrun{
#'  \Sexpr[results=rd,stage=render]{RcppOctave::o_help(rand, format='rd')}
#' }
#' 
#' # to see the included Rd code
#' o_help(rand, format=TRUE)
#' o_help(rand, format='HTML')
#' o_help(rand, format='latex')
#' 
o_help <- function(NAME, character.only = FALSE, show = interactive(), format = c('plain', 'rd', 'txt', 'latex', 'HTML')){
	
	# substitute NAME
	if( !character.only )
		NAME <- as.character(substitute(NAME))
	
	# check argument length
	if( length(NAME) != 1 )
		stop("Argument `NAME` must be a single symbol or character string")
		
    if( isTRUE(format) ) format <- 'rd'
	format <- match.arg(format)
	
	# get the help page from Octave
    hlp <- ''
    compat <- .isPlatformCompatible(details = TRUE)
    if( compat$ok ) hlp <- .Call('oct_help', NAME, PACKAGE='RcppOctave')
    else{
        return(compat$msg)
        # call on platform i386
        libpath <- dirname(path.package('RcppOctave'))
        cmd <- sprintf("Rscript --arch i386 -e \"suppressWarnings(suppressMessages(library(RcppOctave, lib = '%s'))); cat(o_help('%s', TRUE, FALSE, '%s'))\""
                        , libpath
                        , NAME, format)
        return(shell(cmd, intern = TRUE))
    }
    
	#print(hlp)
	if( format != 'plain' ){
		# generate \Sexpr commands for each 
#		sexpr <- paste("RcppOctave::o_help(", NAME, ", show=FALSE)", sep='')
#		rdres <- paste(
#		"\\if{html}{\\Sexpr[results=rd,stage=render]{paste(\"\\\\\\\\out{<pre>\",", sexpr, ",'</pre>}', sep='')}}\n"
#		, "\\if{text}{\\Sexpr[results=rd,stage=render]{paste(\"\\\\\\\\out{\",", sexpr, ", '}', sep='')}}\n"
#		, "\\if{latex}{\\Sexpr[results=rd,stage=render]{paste(\"\\\\\\\\out{\\\\\\\\begin{verbatim}\",", sexpr, ", \"\\\\\\\\end{verbatim}}\", sep='')}}\n"
#		, sep='')
		rdres <- paste(
		"\\if{html}{\\out{<pre>", hlp, "</pre>}}\n"
		, "\\if{text}{\\preformatted{", hlp, "}}\n"
		, "\\if{latex}{\\out{\\begin{verbatim}", hlp, "\\end{verbatim}}}\n"
		, sep='')

		rdo <- c('txt', 'latex', 'HTML')
		if( format=='rd' ) return(rdres)
		else{
			rdfun <- paste('Rd2', format, sep='')
			rdfile <- tempfile("OctaveDoc", fileext=".Rd")
			# write rd code
			write(rdres, file=rdfile)
			# convert in place to required format
			formfile <- tempfile("OctaveDoc", fileext=paste(".", format, sep=''))
			on.exit( unlink(formfile) )
			do.call(rdfun, list(Rd=rdfile, out=formfile, fragment=TRUE))
			formres <- paste(readLines(formfile), collapse="\n")
			if( show ) file.show(formfile, title=, delete.file=TRUE)
			return(invisible(formres))
		}
			
	}
		
	if( show ){
		hlpfile <- tempfile("OctaveDoc", fileext=".txt")
		title <- paste(NAME, "\t\t- Octave Documentation -\t\tv", o_version(), sep='')
		cat(paste(title, "\n\n", hlp, sep=""), file=hlpfile)
		file.show(hlpfile, title=, delete.file=TRUE)
		invisible(hlp)
	}else hlp
}

#' \code{o_doc} displays documentation for the function FUNCTION_NAME directly 
#' from an on-line version of the printed manual, using the GNU Info browser.
#' Type `q` to quit the browser.
#' 
#' @param FUNCTION_NAME the name of the function from which to show the 
#' documentation. See the relevant \emph{Octave Documentation} section below.
#' 
#' @templateVar name doc
#' @template OctaveDoc
#' 
#' @rdname o_help
#' @export
#' @examples
#' 
#' try({ # may throw an error if Octave documentation is not installed
#' 
#' o_doc(text)
#' # or equivalently
#' o_doc('text')
#' 
#' })
#' 
o_doc <- function(FUNCTION_NAME){
	# substitute FUNCTION_NAME
	FUNCTION_NAME <- as.character(substitute(FUNCTION_NAME))
	# call Octave function `doc`
	invisible(.CallOctave('doc', FUNCTION_NAME))
}



#' Octave Identity Function
#' 
#' This function calls the Octave function provided by the module shipped with 
#' RcppOctave. It Returns its arguments unchanged, and is mainly used to test 
#' and check the effect of object conversions between R and Octave. 
#' 
#' @param ... any R object supported by RcppOctave.
#' @return its argument -- list -- after its conversion from R to Octave and 
#' from Octave to R.
#' 
#' @export
#' @examples
#' 
#' o_identity(1L)
#' o_identity(1:10)
#' o_identity(matrix(1:10, 2,5))
#' 
#' o_identity(1)
#' o_identity(runif(10))
#' o_identity(matrix(runif(10), 2,5))
#' 
o_identity <- function(...){
	
	dots <- list(...)
	if( length(dots) > 1 )
		sapply(dots, function(x) .CallOctave('identity', x), simplify=FALSE)
	else
		.CallOctave('identity', dots[[1]])
}

#' \code{o_inpath} tells if a directory or files are in Octave path.
#' 
#' @rdname o_addpath
#' @export
#' @examples
#' 
#' o_addpath(tempdir())
#' o_inpath(tempdir())
#' o_inpath(tempfile())
#'  
o_inpath <- function(...){
	p <- RcppOctave::.CallOctave('path')
    sep <- RcppOctave::.CallOctave('pathsep')
	p <- strsplit(p, sep)[[1]]
	f <- file.path(...)
	sapply(f, function(x){ any(sapply(file.path(p, x), file.exists)) | is.element(x, p)}) 
}

##' Installing An Octave Package
##' 
##' @param x path to a source package
##' @param OPTIONS installation options
##' 
##' @templateVar name pkg
##' @template OctaveDoc
##' 
##' @examples 
##' \dontrun{
##' o_install()
##' }
#o_install <- function(x, OPTIONS=NULL){
#	if( !length(OPTIONS) )
#		.CallOctave('pkg', 'install', x, argout=0)
#	else
#		.CallOctave('pkg', 'install', OPTIONS, x, argout=0)
#	invisible()
#}
