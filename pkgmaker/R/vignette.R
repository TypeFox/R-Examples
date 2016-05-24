# Vignette generation related functions
# 
# Author: Renaud Gaujoux
# Creation: 25 Apr 2012
###############################################################################

#' @include packages.R
NULL

rnw_message <- function(...) message("# ", ...)

#' Identifying Sweave Run
#' 
#' Tells if the current code is being executed within a Sweave document.
#' 
#' @return \code{TRUE} or \code{FALSE}
#' @export
#' 
#' @examples
#' 
#' # Not in a Sweave document 
#' inSweave()
#' 
#' # Within a Sweave document
#' 
inSweave <- function(){
#	in.sweave <- FALSE
	if ((n.parents <- length(sys.parents())) >= 3) {
		for (i in seq_len(n.parents) - 1) {
			if ("chunkopts" %in% ls(envir = sys.frame(i))) {
				chunkopts = get("chunkopts", envir = sys.frame(i))
				if (all(c("prefix.string", "label") %in% names(chunkopts))) {
#					in.sweave <- TRUE
					return(TRUE)
					break
				}
			}
		}
	}
	FALSE
}

#' Generate a Fake Vignette
#' 
#' @param src original Sweave file
#' @param out output file
#' @param PACKAGE package name where to look the source vignette
#' 
#' @export
makeFakeVignette <- function(src, out, PACKAGE=NULL){
	
	# interpret template within the package directory
	if( !is.null(PACKAGE) ){
		src <- str_c(, src)
	}
    if( identical(normalizePath(dirname(src)), normalizePath(dirname(out))) ){
        cat("# NOTE: skipped fake vignette [source in root directory]\n")
        return(invisible())
    }
	# read in template file
	l <- readLines(src)
	# extract %\Vignette commands
	vign <- l[grep("^%\\s*\\\\Vignette", l)]
	# write output file
	cat(c("\\documentclass[10pt]{article}"
		, vign
		, "\\usepackage{url}\n\\usepackage[colorlinks]{hyperref}\n\n\\begin{document}\n\\end{document}")
		, file=out, sep="\n");

}

#' LaTeX Utilities for Vignettes
#' 
#' \code{latex_preamble} outputs/returns command definition LaTeX commands to 
#' be put in the preamble of vignettes.
#' 
#' Argument \code{PACKAGE} is not required for \code{latex_preamble}, but must 
#' be correctly specified to ensure \code{biblatex=TRUE} generates the correct
#' bibliography command.
#'  
#' @param R logical that indicate if general R commands should be added 
#' (e.g. package names, inline R code format commands) 
#' @param CRAN logical that indicate if general CRAN commands should be added
#' (e.g. CRAN package citations) 
#' @param Bioconductor logical that indicate if general Bioconductor commands 
#' should be added (e.g. Bioc package citations) 
#' @param GEO logical that indicate if general GEOmnibus commands should be added
#' (e.g. urls to GEO datasets) 
#' @param ArrayExpress logical that indicate if general ArrayExpress commands 
#' should be added (e.g. urls to ArrayExpress datasets)
#' @param biblatex logical that indicates if a \code{\\bibliography} command 
#' should be added to include references from the package's REFERENCES.bib file. 
#' 
#' @param only a logical that indicates if the only the commands whose 
#' dedicated argument is not missing should be considered.
#' @param file connection where to print. If \code{NULL} the result is returned
#' silently.
#' 
#' @import stringr
#' @export
#' @rdname latex
#' @examples
#' 
#' latex_preamble()
#' latex_preamble(R=TRUE, only=TRUE)
#' latex_preamble(R=FALSE, CRAN=FALSE, GEO=FALSE)
#' latex_preamble(GEO=TRUE, only=TRUE)
#' 
latex_preamble <- function(PACKAGE, R=TRUE, CRAN=TRUE, Bioconductor=TRUE
							, GEO=TRUE, ArrayExpress=TRUE, biblatex=FALSE, only=FALSE, file=''){
	
	cmd <- "%%%% PKGMAKER COMMANDS %%%%%%
\\usepackage{xspace}
"
	
	inc <- function(arg){
		e <- parent.frame()
		(!only || eval(substitute(hasArg(arg), list(arg=substitute(arg))), e)) && arg
	}
	
	if( inc(R) ){
		cmd <- c(cmd, 
"% R
\\let\\proglang=\\textit
\\let\\code=\\texttt 
\\newcommand{\\Rcode}{\\code}
\\newcommand{\\pkgname}[1]{\\textit{#1}\\xspace}
\\newcommand{\\Rpkg}[1]{\\pkgname{#1} package\\xspace}
\\newcommand{\\citepkg}[1]{\\cite{#1}}
")
}

	if( inc(CRAN) ){
		cmd <- c(cmd,
"% CRAN
\\newcommand{\\CRANurl}[1]{\\url{http://cran.r-project.org/package=#1}}
%% CRANpkg
\\makeatletter
\\def\\CRANpkg{\\@ifstar\\@CRANpkg\\@@CRANpkg}
\\def\\@CRANpkg#1{\\href{http://cran.r-project.org/package=#1}{\\pkgname{#1}}\\footnote{\\CRANurl{#1}}}
\\def\\@@CRANpkg#1{\\href{http://cran.r-project.org/package=#1}{\\pkgname{#1}} package\\footnote{\\CRANurl{#1}}}
\\makeatother
%% citeCRANpkg
\\makeatletter
\\def\\citeCRANpkg{\\@ifstar\\@citeCRANpkg\\@@citeCRANpkg}
\\def\\@citeCRANpkg#1{\\CRANpkg{#1}\\cite*{Rpackage:#1}}
\\def\\@@citeCRANpkg#1{\\CRANpkg{#1}~\\cite{Rpackage:#1}}
\\makeatother
\\newcommand{\\CRANnmf}{\\href{http://cran.r-project.org/package=NMF}{CRAN}}
\\newcommand{\\CRANnmfURL}{\\url{http://cran.r-project.org/package=NMF}}
")
}

	if( inc(Bioconductor) ){
		cmd <- c(cmd,
"% Bioconductor
\\newcommand{\\BioCurl}[1]{\\url{http://www.bioconductor.org/packages/release/bioc/html/#1.html}}
\\newcommand{\\BioCpkg}[1]{\\href{http://www.bioconductor.org/packages/release/bioc/html/#1.html}{\\pkgname{#1}} package\\footnote{\\BioCurl{#1}}}
\\newcommand{\\citeBioCpkg}[1]{\\BioCpkg{#1}~\\cite{Rpackage:#1}}
% Bioconductor annotation
\\newcommand{\\BioCAnnurl}[1]{\\url{http://www.bioconductor.org/packages/release/data/annotation/html/#1.html}}
\\newcommand{\\BioCAnnpkg}[1]{\\href{http://www.bioconductor.org/packages/release/data/annotation/html/#1.html}{\\Rcode{#1}} annotation package\\footnote{\\BioCAnnurl{#1}}}
\\newcommand{\\citeBioCAnnpkg}[1]{\\BioCAnnpkg{#1}~\\cite{Rpackage:#1}}
")
}

	if( inc(GEO) ){
		cmd <- c(cmd, 
"% GEO
\\newcommand{\\GEOurl}[1]{\\href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=#1}{#1}\\xspace}
\\newcommand{\\GEOhref}[1]{\\GEOurl{#1}\\footnote{\\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=#1}}}
")
	}

	if( inc(ArrayExpress) ) cmd <- c(cmd,
"% ArrayExpress
\\newcommand{\\ArrayExpressurl}[1]{\\href{http://www.ebi.ac.uk/arrayexpress/experiments/#1}{#1}\\xspace}
\\newcommand{\\ArrayExpresshref}[1]{\\ArrayExpressurl{#1}\\footnote{\\url{http://www.ebi.ac.uk/arrayexpress/experiments/#1}}}
")

	if( biblatex ){
		if( missing(PACKAGE) )
			stop("Argument `PACKAGE` is required when specifying `biblatex=TRUE`.")
		cmd <- c(cmd, latex_bibliography(PACKAGE, file=NULL))
	} 

	# output or return commands
	cmd <- c(cmd, "%%%% END: PKGMAKER COMMANDS %%%%%%\n")
	cmd <- str_c(cmd, collapse="\n")
	if( !is.null(file) ) cat(cmd, file, sep='')
	else cmd
	
}

#' \code{latex_bibliography} prints or return a LaTeX command that includes a 
#' package bibliography file if it exists.
#' 
#' @param PACKAGE package name
#' 
#' @export
#' @rdname latex
#' 
latex_bibliography <- function(PACKAGE, file=''){
	
    rpkg.bib <- "%\\bibliography{Rpackages}\n"
    cmd <- rpkg.bib
	# get REFERENCES.bib file
	reffile <- packageReferenceFile(PACKAGE=PACKAGE)
    if( is.file(reffile) ){
        cmd <- paste0(cmd, "\\bibliography{", gsub("\\.bib$", "", reffile), "}\n")
    }
	
    # add post-processing knit hook
    library(knitr)
    knit_hooks$set(document = function(x){
        # write bibfile if necessary
        if( length(pkgs <- parsePackageCitation(x)) ){
            # write bibfile
            write.pkgbib(gsub("^Rpackage:", '', pkgs), file='Rpackages.bib', prefix='Rpackage:')
            # uncomment inclusion line
            x <- gsub("%\\bibliography{Rpackages}", "\\bibliography{Rpackages}", x, fixed = TRUE)
        }
        x
    })

    if( !is.null(file) ) cat(cmd, file=file)
    else cmd
}

is.rnw <- function(x){
	is(x, 'rnw')
}


runVignette <- function(x, ...){
    # flag the vignette as being locally made
    opts <- options(R_RUNNING_MANUAL_VIGNETTE=TRUE)
    on.exit( options(opts) )
    # run
	UseMethod('runVignette')
}

# tells if a vignette is locally made
#' Identifies Manually Run Vignettes
#' 
#' \code{isManualVignette} tells if a vignette is being run through the function \code{runVignette} 
#' of \pkg{pkgmker}, allowing disabling behaviours not allowed in package vignettes that are
#' checked vi \code{R CMD check}. 
#' 
#' @rdname vignette
#' @export
isManualVignette <- function(){
    isTRUE(getOption('R_RUNNING_MANUAL_VIGNETTE'))
}

runVignette.default <- function(x, file=NULL, ...){
	stop("Vignette compiler '", class(x), "' is not supported")
} 

## #' @param fig.path specification for the figure path (used in knitr vignettes only). 
## #' If \code{TRUE} then the figure path is set to \code{'./figure/<basename>/'}.
## #' @param cache.path specification for the cache path.
## #' If \code{TRUE} then the figure path is set to \code{'./cache/<basename>/'}.
#' @S3method runVignette rnw_knitr 
runVignette.rnw_knitr <- function(x, file=NULL, ..., fig.path=TRUE, cache.path=TRUE){
	library(knitr)
	# expand path to cache to fix issue in knitr
	bname <- sub("\\..{3}$", '', basename(x$file))	
	# add suffix for windows
	if( .Platform$OS.type == 'windows' ){ 
		bname <- paste(bname, '-win', sep='')
	}
	# cache.path
	if( !isFALSE(cache.path) ){
		if( isTRUE(cache.path) ){
			cache.path <- file.path(getwd(), 'cache', bname, '/')
		}
		opts_chunk$set(cache.path=cache.path)	
	}
	# fig.path
	if( !isFALSE(fig.path) ){
		if( isTRUE(fig.path) ){
			fig.path <- file.path(getwd(), 'figure', str_c(bname,'-'))
		}
		opts_chunk$set(fig.path=fig.path)	
	}
	
	# set other options
	opts_chunk$set(...)
	
	# run knitr
	e <- new.env(parent = .GlobalEnv)
	if( FALSE && (is.null(file) || file_extension(file) %in% c('tex', 'pdf')) ){
		ofile <- if( file_extension(file) == 'pdf' ) file else NULL 
		knit2pdf(x$file, ofile, envir=e)
		if( is.null(file) ){
			# remove pdf file
			unlink(file.path(getwd(), basename(file_extension(x$file, 'pdf'))))
		} else if( file_extension(file) == 'tex' ){
			# move tex file
			file.rename(file_extension(file, 'tex'), file)
		}
	}else knit(x$file, file, envir=e)
}

#' @S3method runVignette rnw_sweave 
runVignette.rnw_sweave <- function(x, file=NULL, ...){
	res <- Sweave(x$file, driver=x$driver, ...)
	# move output file
	if( !is.null(file) ){
		file.rename(res, file)
	}
	res
}

#' Utilities for Vignettes
#' 
#' \code{rnw} provides a unified interface to run vignettes that detects
#' the type of vignette (Sweave or \code{\link[knitr]{knitr}}), and which Sweave driver 
#' to use (either automatically or from an embedded command \code{\\VignetteDriver} 
#' command).
#' 
#' @param x vignette source file specification as a path or a \code{rnw} object.
#' @param file output file
#' @param ... extra arguments passed to \code{as.rnw} that can be used to force 
#' certain building parameters.
#' @param raw a logical that indicates if the raw result for the compilation 
#' should be returned, instead of the result file path.
#'   
#' @rdname vignette
#' @export
rnw <- function(x, file=NULL, ..., raw=FALSE){
	
#	library(methods)
	# load rnw file
	x <- as.rnw(x, ...)	
	
	# setup restoration of options and graphical parameters
	opts <- options()
	gpar <- par(no.readonly=TRUE)
	on.exit( {options(opts); par(gpar)})
	
	# copy package cleveref from pkgmaker installation
	if( 'cleveref' %in% x$latexPackages ){
		clv <- packagePath('cleveref.sty', package='pkgmaker')
		message("# Copying package 'cleveref.sty' from ", dirname(clv)," ... ", appendLF=FALSE)
		wd <- if( !is.null(file) ) dirname(file) else getwd()
		file.copy(clv, wd)
		if( file.exists(file.path(wd, basename(clv))) )	message('OK') else message('ERROR')
	}
	
	# run vignette
	res <- runVignette(x, file=file, ...)
	
	# Package citations
	if( !is.null(keys <- x$cite) ){
		message("# Writing package bibtex file [", length(keys)," key(s)] ... ", appendLF=FALSE)
		write.pkgbib(keys, file='Rpackages.bib', prefix='Rpackage:', verbose=FALSE)
		message('OK')
	}
	#

	# return raw result if required
	if( raw ) return(res)

	# check for a wrapper main file
	if( !is.null(x$wrapper) ){
		res <- x$wrapper
	}
	
	invisible(res)
}

checkFile <- function(x, msg="file '%s' does not exist."){
	if( !is.file(x) ) stop(sprintf(msg, x))
	TRUE
}

checkRnwFile <- function(x){
	if( is.rnw(x) ) x <- x$file
	checkFile(x, msg="Vignette file '%s' does not exist.")
}

#' \code{as.rnw} creates a S3 \code{rnw} object that contains information
#' about a vignette, e.g., source filename, driver, fixed included files, etc..
#' 
#' @param load logical to indicate if all the object's properties should loaded, 
#' which is done by parsing the file and look up for specific tags. 
#' 
#' @rdname vignette
#' @export
as.rnw <- function(x, ..., load = TRUE){
	
	if( is.rnw(x) ) return(x)
	
	checkRnwFile(x)
	# initialise 'rnw' object
	obj <- list()
	class(obj) <- 'rnw'
	
	# store source full path
	obj$file <- normalizePath(x)
	obj$line <- NA
	if( !load ) return(obj)
	
	# detect compiler
	obj$compiler <- rnwCompiler(obj) %||% 'Sweave'
	cl <- if( obj$compiler == 'knitr' ) 'knitr' else 'sweave'
	class(obj) <- c(paste('rnw_', cl, sep=''), class(obj))
	
	# detect driver
	obj$driver <- rnwDriver(obj) %||% RweaveLatex()
	# detect wrapper
	obj$wrapper <- rnwWrapper(obj)
	# detect fixed included images
	obj$includes <- rnwIncludes(obj)
	# detect latex packages
	obj$latexPackages <- rnwLatexPackages(obj)
	# detect children vignettes
	obj$children <- rnwChildren(obj)
	# detect package citations
	obj$cite <- rnwCite(obj)
	
	# override with passed extra arguments
	if( nargs() > 1L ){
		dots <- list(...)
		obj[names(dots)] <- dots
	}
	
	# return object
	obj
} 

rnwObject <- function(...) as.rnw(..., load=FALSE)

rnwParser <- function(tag, name=tolower(tag), trim=TRUE, commented=FALSE, options=FALSE, first=FALSE){
	
	function(x, verbose=TRUE){
		x <- rnwObject(x)
		# read all lines in
		l <- readLines(x$file)
		
		# identify driver
		dr <- str_match(l, str_c("^\\s*"
								, if( commented ) '%'
								,"\\s*\\\\", tag
								, if( options ) "(\\[[^]]*\\])?"
								, "\\{([^}]*)\\}"))
		w <- which(!is.na(dr[,1L]))
		if( length(w) > 0L ){
			if( first ) w <- w[1L]
			s <- dr[w, if( options ) 3L else 2L]
			# trim if necessary
			if( trim ) s <- str_trim(s)
			if( verbose ) rnw_message("Detected ", name, ": "
								,paste("'", s, "'", sep='', collapse=', '))
			s
		}
	}
}

rnwVignetteParser <- function(tag, ...){
	rnwParser(str_c('Vignette',tag), name=tolower(tag), ..., commented=TRUE, first=TRUE)
}

rnwLatexPackages <- rnwParser('usepackage', name='LaTeX package(s)', options=TRUE)

#' \code{rnwCompiler} tries to detect the vignette compiler to use on a vignette
#' source file, e.g., \code{\link{Sweave}} or \code{\link[knitr]{knitr}}.
#' 
#' @param verbose logical that toggles verbosity
#' 
#' @rdname vignette
#' @export
rnwCompiler <- rnwVignetteParser('Compiler')

#' \code{rnwWrapper} tries to detect the type of vignette and if it is meant 
#' to be wrapped into another main file.
#' 
#' @rdname vignette
#' @export
rnwWrapper <- rnwVignetteParser('Wrapper')

#' \code{rnwDriver} tries to detect Sweave driver to use on a vignette source 
#' file, e.g., \code{SweaveCache}, \code{highlight}, etc..
#' 
#' @rdname vignette
#' @export
rnwDriver <- function(x){
	
	parse_driver <- rnwVignetteParser('Driver', trim=FALSE)
	if( !is.null(s <- parse_driver(x)) ){
		# eval text
		eval(parse(text=s))
	}
	
}  

#' \code{rnwIncludes} detects fixed includes, e.g., image or pdf files, that are 
#' required to build the final document.  
#' 
#' @rdname vignette
#' @export
rnwIncludes <- function(x){
	
	x <- rnwObject(x)
	
	# read all lines in
	l <- readLines(x$file)
	
	# identify driver
	dr <- suppressWarnings(str_match(l, "^\\s*\\\\((include)|(includegraphics)|(input))\\{([^}]*)\\}"))
	w <- which(!is.na(dr[,1L]))
	rnw_message("Detected includes: ", appendLF=FALSE)
	if( length(w) > 0L ){
		inc <- str_trim(dr[w,6L])
		message(str_out(inc))
		inc
	}else
		message("NONE")
	
}

#' \code{rnwChildren} detects included vignette documents and return them as a 
#' list of vignette objects.  
#'  
#' @rdname vignette
#' @export
rnwChildren <- function(x){
	
	x <- rnwObject(x)
	
	# read all lines in
	l <- readLines(x$file)
	
	# identify driver
	dr <- str_match(l, "^\\s*\\\\SweaveInput\\{([^}]*)\\}")
	w <- which(!is.na(dr[,1L]))
	if( length(w) > 0L ){
		inc <- dr[w,2L]
		rnw_message("Detected children: ", str_out(inc, Inf))
		owd <- setwd(dirname(x$file))
		on.exit( setwd(owd) )
		mapply(as.rnw, inc, line=w, SIMPLIFY=FALSE)
	}
	
}  

#' Formatting Package Citations in Sweave/knitr Documents
#' 
#' @param x output document, as a single string.
#' @export
parsePackageCitation <- function(x){
    
    if( length(x) > 1L ) x <- paste(x, collapse = "\n")
    
    .parse <- function(x, pattern, idx){
		dr <- str_match_all(x, pattern)
		dr <- dr[sapply(dr, length)>0L]
		unlist(lapply(dr, '[', , idx))
	}
	
	# extract package citations: \citeCRANpkg - like
    x <- gsub(".*[^%]* *\\\\begin\\{document\\}(.*)", "\\1", x)
    cite <- .parse(x, "\\\\cite((CRAN)|(BioC)|(BioCAnn))?pkg[*]?\\{([^}]*)\\}", 6L)
    # \cite{Rpackage:pkgname, ...} - like
	cite2 <- .parse(x, "\\\\cite[^{ ]*\\{([^}]*)\\}", 2L)
    if( length(cite2) ){
 		cite2 <- .parse(cite2, '.*Rpackage:([^,}]+).*', 2L)
		cite <- c(cite, cite2)
	}
	# remove Rpackage prefix
	if( length(cite) ){
        cite <- unlist(strsplit(cite, ","))
        cite <- gsub('^Rpackage:', '', cite)
	}
	
    inc <- character()
	if( length(cite) > 0L ){
		inc <- unique(str_trim(unlist(strsplit(cite, ","))))
    }
	inc  
}

#' \code{bibcite} provides an inline package citation functionnality. 
#' Technically it adds a given Bibtex key to a cache that is used at the end of the 
#' document processing to generate a .bib file with all citation keys.  
#' 
#' @param key citation Bibtex key(s) as a character vector
#' @param cache specifies what to do with the previsouly chached keys.
#' If \code{TRUE}, then \code{key} is added to the cache. 
#' If \code{NULL}, then all previously cached keys are deleted, before .
#' If a character string, then it specifies the path to a Bibtex file that is loaded 
#' to initialise the cache.
#' @param ... extra arguments passed to \code{\link[bibtex]{read.bib}}.
#' @keywords internal
cite_pkg <- local({
    .keys <- character()
    function(key, cache = NA, ...){
        # return whole cache
        if( !nargs() ) return(.keys)
        # reset cache
        if( is.null(cache) ) .keys <- character()
        else if( isString(cache) ) .keys <- read.bib(file = cache, ...)
        if( !missing(key) ){
            cat(key)
            .keys <<- c(.keys, key)
        }
    }
})

rnwCite <- function(x){
	
	x <- rnwObject(x)
	
	# read all lines in
	l <- readLines(x$file)

	.parse <- function(x, pattern, idx){
		dr <- str_match_all(x, pattern)
		dr <- dr[sapply(dr, length)>0L]
		unlist(lapply(dr, '[', , idx))
	}
	
	# extract package citations: \citeCRANpkg - like
	cite <- .parse(l, "\\\\cite((CRAN)|(BioC)|(BioCAnn))?pkg[*]?\\{([^}]*)\\}", 6L)
	# \cite{Rpackage:pkgname, ...} - like
	cite2 <- .parse(l, "\\\\cite[^{ ]*\\{([^}]*)\\}", 2L)
	if( length(cite2) ){
 		cite2 <- .parse(cite2, '.*Rpackage:([^,}]+).*', 2L)
		cite <- c(cite, cite2)
	}
	# remove Rpackage prefix
	if( length(cite) ){
    cite <- unlist(strsplit(cite, ","))
	  cite <- gsub('^Rpackage:', '', cite)
	}
	
	rnw_message("Detected package citation(s): ", appendLF=FALSE)
	if( length(cite) > 0L ){
		inc <- unique(str_trim(unlist(strsplit(cite, ","))))
		message(str_out(inc), ' [', length(inc), ']')
		inc
	}else
		message("NONE")
}

# substitute a makefile template variable 
subMakeVar <- function(mvar, value, text){
	sub(str_c('#%', mvar, '%#'), value, text, fixed=TRUE)
}
# define a makefile template variable
defMakeVar <- function(var, value, ..., mvar=var){
	subMakeVar(mvar, str_c(var, '=', value), ...)
}

quick_install <- function(path, ..., lib.loc){
	
	if( !is.dir(lib.loc) ){
		stop("Installation directory '", lib.loc, "' does not exist.")
	}
	
	olib <- .libPaths()
	.libPaths(lib.loc)
	on.exit( .libPaths(olib) )
	pkg <- devtools::install(path, ...)
	
}


vignetteCheckMode <- checkMode_function('_R_CHECK_BUILDING_VIGNETTES_')

#' \code{vignetteMakefile} returns the path to a generic makefile used to make 
#' vignettes.
#' 
#' @param package package name.
#' If \code{NULL}, a DESRIPTION file is looked for one directory up: this 
#' meant to work when building a vignette directly from a package's 
#' \code{'vignettes'} sub-directory. 
#' @param skip Vignette files to skip (basename).  
#' @param print logical that specifies if the path should be printed or
#' only returned. 
#' @param template template Makefile to use.
#' The default is to use the file \dQuote{vignette.mk} shipped with the package
#' \pkg{pkgmaker} and can be found in its install root directory.
#' @param temp logical that indicates if the generated makefile should using a 
#' temporary filename (\code{TRUE}), or simply named \dQuote{vignette.mk}
#' @param checkMode logical that indicates if the vignettes should be generated as in a 
#' CRAN check (\code{TRUE}) or in development mode, in which case \code{pdflatex}, \code{bibtex}, 
#' and, optionally, \code{qpdf} are required.
#' @param user character vector containing usernames that enforce \code{checkMode=TRUE}, 
#' if the function is called from within their session.
#' @param tests logical that enables the compilation of a vignette that gathers all unit 
#' test results.
#' Note that this means that all unit tests are run before generating the vignette.
#' However, unit tests are not (re)-run at this stage when the vignettes are built 
#' when checking the package with \code{R CMD check}.
#' 
#' @rdname vignette
#' @export
vignetteMakefile <- function(package=NULL, skip=NULL, print=TRUE, template=NULL, temp=FALSE
                             , checkMode = isCHECK() || vignetteCheckMode()
                             , user = NULL, tests=TRUE){
	
#	library(methods)
	## create makefile from template
	# load template makefile
	if( is.null(template) )
		template <- packagePath('vignette.mk', package='pkgmaker')	
	l <- paste(readLines(template), collapse="\n")
  
	# R_BIN
	l <- subMakeVar('R_BIN', R.home('bin'), l)
  #
  
    if( checkMode ){
        oldCM <- vignetteCheckMode(TRUE)
        on.exit( vignetteCheckMode(oldCM) )
    }
  # Check user: LOCAL_MODE if in declared user
	localMode <- !checkMode
	cuser <- Sys.info()["user"]
	l <- subMakeVar('VIGNETTE_USER', cuser, l)
  maintainers <- '-'
  if( !is.null(user) ){
		maintainers <- str_c(user, collapse=', ')
		if( cuser %in% user ){
			localMode <- TRUE
		}
	}
	l <- subMakeVar('VIGNETTE_MAINTAINERS', maintainers, l)
	
  # define variable LOCAL_MODE
	if( localMode ){
	  l <- defMakeVar('LOCAL_MODE', cuser, l)
  }
  
	# Package name
    pkg_dir <-  dirname(getwd())
    loc_package <- if( is.file(df <- file.path(pkg_dir, 'DESCRIPTION')) ){
        d <- try(read.dcf(df), silent=TRUE)
        d <- as.list(as.data.frame(d, stringsAsFactors=FALSE))
        d$Package
    }
    if( !is.null(loc_package) && (is.null(package) || identical(loc_package, package)) ) package <- loc_package
    else if( !identical(loc_package, package) && length(pkg_dir <- find.package(package, quiet=TRUE)) ){
		d <- packageDescription(package)
	}else{
		stop("Could not load DESCRIPTION file for package '", package, "'.")		
	}
	l <- defMakeVar('MAKE_R_PACKAGE', package, l)
    l <- subMakeVar('R_PACKAGE_DESCRIPTION', pkg_dir, l)
  # R_LIBS: add package's dev lib if necessary 
  Rlibs <- NULL
  if( localMode && is.dir(devlib <- file.path(getwd(), '..', '..', 'lib')) ){
	  Rlibs <- devlib
  }
  Rlibs <- paste(c(Rlibs, "$(TMP_INSTALL_DIR)", "$(R_LIBS)"), collapse=.Platform$path.sep)
	l <- subMakeVar('R_LIBS_DEV', Rlibs, l)
	# TMP_INSTALL_DIR: temporary install directory
	l <- subMakeVar('TMP_INSTALL_DIR', file.path(dirname(tempdir()), basename(tempfile('Rpkglib_'))), l)
	
	# Vignettes files:
    # - look into src/ for real vignettes
	# - check presence of a test directory ../tests/
	# - check current directory for non fake vignettes
	rnwFiles <- NULL
	# src
	if( is.dir('src') ) rnwFiles <- list.files('src', pattern="\\.Rnw$")
	# unit tests
	if( tests && is.dir('../tests') ) rnwFiles <- c(rnwFiles, str_c(package, '-unitTests.Rnw'))
	# non-fake vignettes
    rnwFiles <- c(rnwFiles, list.files('.', pattern="\\.Rnw$"))
	# substitute in makefile
	rnwFiles <- unique(rnwFiles)
	if( !is.null(skip) )
		rnwFiles <- setdiff(rnwFiles, skip)
	l <- subMakeVar('RNW_SRCS', paste(rnwFiles, collapse=' '), l)
	# reset pdf objects in local mode to point to ../inst/doc
	noBuildVignettes <- if( !is.null(d$BuildVignettes) ) tolower(d$BuildVignettes)=='no' else FALSE
	if( localMode && noBuildVignettes ){
        l <- defMakeVar('INST_TARGET', 1, l)
    	l <- defMakeVar('PDF_OBJS'
						, paste(file.path('../inst/doc', sub("\\.Rnw$", ".pdf", rnwFiles)), collapse=' ')
						, l)
	}
    l <- defMakeVar('PDF_OBJS'
            , paste(file.path('../inst/doc', sub("\\.Rnw$", ".pdf", rnwFiles)), collapse=' ')
            , l)
    
	# create makefile
	mk <- if( temp ) tempfile('vignette_', tmpdir='.', fileext='.mk') else 'vignette.mk'
	cat(l, file=mk)
	if ( print ){
		cat(mk)
	}
	invisible(l)
}

#' Compact PDF at Best
#' 
#' Compact PDFs using either \code{gs_quality='none'} or \code{'ebook'}, 
#' depending on which compacts best (as per CRAN check criteria).
#' 
#' @inheritParams tools::compactPDF
#' 
#' @rdname vignette
#' @export
compactVignettes <- function(paths, ...){
	
	td <- tempfile(basename(paths))
	file.copy(paths, td)
	res <- tools::compactPDF(td, gs_quality = "none", ...) # use qpdf
	diff_none <- format(res, diff = 1e5)
	res <- tools::compactPDF(td, gs_quality = "ebook", ...)
	diff_ebook <- format(res, diff = 2.5e5) # 250 KB for now
	
	if( length(diff_ebook) ){
		tools::compactPDF(paths, gs_quality = "ebook", ...)
		invisible('ebook')
	}else{
		tools::compactPDF(paths, gs_quality = "none", ...)
		invisible('none')
	}
	
}

