# Package related functions
# 
# Author: Renaud Gaujoux
# Creation: 29 Jun 2012
###############################################################################

path.protect <- function(...){
  f <- file.path(...)
  if( .Platform$OS.type == 'windows' ){
    f <- gsub("\\\\", "/", f)
  }
  paste('"', f, '"', sep='')
}

#' Quick Installation of a Source Package
#' 
#' Builds and install a minimal version of a package from its 
#' source directory.
#' 
#' @param path path to the package source directory
#' @param destdir installation directory. 
#' If \code{NULL}, the package is installed in the default installation library.
#' If \code{NA}, the package is installed in a temporary directory, whose path is returned
#' as a value.
#' @param vignettes logical that indicates if the vignettes should be 
#' rebuilt and installed.
#' @param force logical that indicates if the package should be installed even if a previous
#' installation exists in the installation library.
#' @param ... extra arguments passed to \code{\link{R.CMD}}
#' @param lib.loc library specification.
#' If \code{TRUE} then the installation directory \code{destdir} is added to the default 
#' library paths.
#' This can be usefull if dependencies are installed in this directory.
#' If \code{NULL}, then the default library path is left unchanged.
#' 
#' @return The path of the library where the package was installed.
#' 
#' @export
#' 
quickinstall <- function(path, destdir=NULL, vignettes=FALSE, force=TRUE, ..., lib.loc=if(!is.null(destdir)) TRUE){
	
	npath <- normalizePath(path)
	pkg <- as.package(path)
  
  # define installation library
	nlib <- if( !is.null(destdir) ) destdir
  else if( is_NA(destdir) ) tempfile("pkglib_")
  
  # normalize path
  if( !is.null(nlib) ){
    # create direcory if needed
    if( !is.dir(nlib) ) dir.create(nlib, recursive=TRUE)
    nlib <- normalizePath(nlib)
    
    if( !is.dir(nlib) ){
      stop("Could not install package '", pkg$package, "': installation directory '", nlib, "' does not exist.")
    }
    
    # add destination directory to default libraries
    if( isTRUE(lib.loc) ) lib.loc <- unique(c(nlib, .libPaths()))
  }
  
  # setup result string
	res <- invisible(if( !is.null(destdir) ) nlib else .libPaths()[1L])
  
  # early exit if the package already exists in the library (and not forcing install)
	message("# Check for previous package installation ... ", appendLF=FALSE)
  if( !is.null(destdir) && is.dir(file.path(nlib, pkg$package)) ){
    if( !force ){
      message("YES (skip)")
      return(res)
    }
    message("YES (replace)")
  }else message("NO")
  
	# add lib path
	ol <- set_libPaths(lib.loc)
	on.exit(set_libPaths(ol), add=TRUE)
	message("Using R Libraries: ", str_out(.libPaths(), Inf))
	
	owd <- setwd(tempdir())
	on.exit( setwd(owd), add=TRUE)
  
	# build
	message("# Building package `", pkg$package, "` in '", getwd(), "'")
	opts <- '--no-manual --no-resave-data '
	if( !vignettes ){
        vflag <- if( testRversion('>= 3.0') ) '--no-build-vignettes ' else '--no-vignettes ' 
        opts <- str_c(opts, vflag)
    }
	R.CMD('build', opts, path.protect(npath), ...)
	spkg <- paste(pkg$package, '_', pkg$version, '.tar.gz', sep='')
	if( !file.exists(spkg) ) stop('Error in building package `', pkg$package,'`')
	# install
	message("# Installing package `", pkg$package, "`"
          , if( !is.null(destdir) ){
            tmp <- if( is_NA(destdir) ) 'temporary '
            str_c("in ", tmp, "'", nlib, "'")
          })
  opts_inst <- ' --no-multiarch --no-demo --with-keep.source '
	if( !vignettes ) opts_inst <- str_c(opts_inst, '--no-docs ')
	R.CMD('INSTALL', if( !is.null(destdir) ) paste('-l', path.protect(nlib)), opts_inst, path.protect(spkg), ...)
  
  # return installation library
	invisible(res)
}

#' Require a Package
#' 
#' Require a package with a custom error message
#' 
#' @param pkg package name as a character string
#' @param ... extra arguments concatenated to for the header of the 
#' error message 
#' 
#' @export
requirePackage <- function(pkg, ...){
	
	if( !require(pkg, character.only=TRUE) ){
		if( nargs() > 1L ) stop(..., " requires package(s) ", str_out(pkg))
		else stop("Could not find required package(s) ", str_out(pkg))
	}
}

parse_deps <- function (string) 
{
	if (is.null(string)) 
		return()
	string <- gsub("\\s*\\(.*?\\)", "", string)
	pieces <- strsplit(string, ",")[[1]]
	pieces <- gsub("^\\s+|\\s+$", "", pieces)
	pieces[pieces != "R"]
}

#' List Package Dependencies
#' 
#' @param x path to package source directory or file.
#' @param all logical that indicates if all dependencies should be returned,
#' or only the required ones.
#' @param as.list logical that indicates if the result should be a list with one element
#' per type of dependency.
#' @param available a matrix of available packages (as returned by \code{\link{available.packages}}), 
#' from which the dependencies are retrieved.
#' This means that there must be a row for the package \code{x}.
#'  
#' @export
#' 
packageDependencies <- function(x, all = TRUE, as.list = FALSE, available = NULL){
    
    if( is.null(available) ) x <- as.package(x, extract = TRUE)
    else{
        p <- available[, 'Package']
        if( !x %in% p ) return(NA)
        x <- available[p == x, , drop = FALSE][1L, ]
        names(x) <- tolower(names(x))
    }
    
	d <- lapply(x[c('depends', 'imports', 'linkingto', 'suggests')], parse_deps)
	d <- unlist(d)
    d <- d[!is.na(d)]
    if( !length(d) ) return()
    names(d) <- gsub("[0-9]+$", "", names(d))
    
    # remove non-required
    if( !all ) d <- d[!names(d) %in% c('suggests')]    
    if( as.list ) d <- split(unname(d), names(d))
    d
}

.biocLite <- function(...){
	# install BiocInstaller if necessary
	if( !require.quiet('BiocInstaller') ){
		message("Installing biocLite")
		source('http://www.bioconductor.org/biocLite.R')
	}
	f <- get('biocLite', 'package:BiocInstaller')
	f(...)
}

#' Installing All Package Dependencies
#' 
#' Install all dependencies from a package source directory or 
#' package source file. 
#' 
#' @param pkg package path or source file
#' @param all logical that indicates if 'Suggests' packages
#' should be installed.
#' @param ... extra arguments passed to \code{\link{install.packages}}.
#' @param dryrun logical that indicates if the packages should be 
#' effectively installed or only shown. 
#' 
#' @export
#' @examples 
#' 
#' try( install.dependencies('Matrix', dryrun=TRUE) )
#' \dontrun{
#' install.dependencies("mypackage_1.0.tar.gz", dryrun=TRUE)
#' }
#' 
install.dependencies <- function (pkg = NULL, all=FALSE, ..., dryrun=FALSE) 
{
	pkg <- as.package(pkg, extract=TRUE)
	deps <- c(parse_deps(pkg$depends)
			, parse_deps(pkg$imports) 
			, parse_deps(pkg$linkingto)
			, if( isTRUE(all) ) parse_deps(pkg$suggests) )
	not.installed <- function(x) length(find.package(x, quiet = TRUE)) == 0
	message("Package dependencies for ", pkg$package, ": ", str_out(deps, Inf))
	deps <- Filter(not.installed, deps)
	if (length(deps) == 0){
		message("Missing: none")
		return(invisible())
	}
	message("Missing: ", str_out(deps, Inf))
	message("Installing ", length(deps), " dependencies for ", pkg$package)
	if( !dryrun ){
		.biocLite(deps, ...)
	}
	invisible(deps)
}

#' Setting Mirrors and Repositories
#' 
#' \code{setBiocMirror} sets all Bioconductor repositories (software, data, 
#' annotation, etc.).
#' so that they are directly available to \code{\link{install.packages}}.
#' It differs from \code{\link{chooseBioCmirror}} in that it effectively enables 
#' the repositories.
#' 
#' @param url or Bioconductor mirror url
#' @param version version number
#' @param unique logical that indicate if duplicated urls or names should be 
#' removed.
#'
#' @rdname mirrors
#' @export 
setBiocMirror <- function(url='http://www.bioconductor.org', version=NULL, unique=TRUE){
	
    #get all bioconductor repos      
    biocRepos <- getBiocRepos(url, version)
	
	repos <- c(biocRepos, getOption('repos'))
	if( unique ){
		nam <- names(repos)
		repos <- repos[!duplicated(repos) & (!duplicated(nam) | nam=='')]
	}
    options(repos=repos)
}

#' \code{getBiocMirror} is a shortcut for \code{getOption('BioC_mirror')}, which 
#' returns the current Bioconductor mirror as used by \code{biocLite}.
#'  
#' @export
#' @rdname mirrors
getBiocMirror <- function(){
	getOption('BioC_mirror')
}
#' \code{getBiocRepos} returns urls to all Bioconductor repositories on a 
#' given mirror.
#' 
#' @export
#' @rdname mirrors
getBiocRepos <- function(url='http://www.bioconductor.org', version=NULL){
	
	if( is.null(url) ){
		url <- getBiocMirror()
		if( is.null(url) )
			stop("No Bioconductor mirror was setup. Use `setBiocMirror`.")
	}
	
	## BioConductor CRAN-style repositories.
	## The software repo (bioc) _must_ be the first element.
	biocParts <- c(
			bioc='bioc'
			, biocData='data/annotation'
			, biocExp='data/experiment'
			, biocExtra='extra'
    )
	
	# define version suffix for bioconductor repo
	if( is.null(version) ){
		assoc <- list(`2`=c(7L, 2L))
		Rv <- as.integer(sub("([0-9]+).*", "\\1", R.version$minor))
		offset <- assoc[[R.version$major]]
	    version <- paste(R.version$major, offset[2L] + Rv - offset[1L], sep='.')
	}
	
	#add version suffix for bioconductor repo
    setNames(paste(url, 'packages', version, biocParts, sep='/'), names(biocParts))
}

#' \code{setCRANMirror} sets the preferred CRAN mirror.
#' 
#' @rdname mirrors
#' @export
setCRANMirror <- function(url=CRAN, unique=TRUE){
	
	repos <- c(CRAN=url, getOption('repos'))
	if( unique ){
		nam <- names(repos)
		repos <- repos[!duplicated(repos) & (!duplicated(nam) | nam=='')]
	}
    options(repos=repos)
}

#' \code{CRAN} simply contains the url of CRAN main mirror 
#' (\url{http://cran.r-project.org}), and aims at simplifying its use, e.g., in 
#' calls to \code{\link{install.packages}}.
#' 
#' @rdname mirrors
#' @export
#' 
#' @examples
#' \dontrun{
#' install.packages('pkgmaker', repos=CRAN)
#' }
CRAN <- 'http://cran.r-project.org'


#' Adding Package Libraries
#' 
#' Prepend/append paths to the library path list, using \code{\link{.libPaths}}.
#' 
#' This function is meant to be more convenient than \code{.libPaths}, which requires 
#' more writing if one wants to:
#' \itemize{
#' \item sequentially add libraries;
#' \item append and not prepend new path(s);
#' \item keep the standard user library in the search path.
#' }
#' 
#' @param ... paths to add to .libPath
#' @param append logical that indicates that the paths should be appended
#' rather than prepended.
#' 
#' @export
#' 
#' @examples
#' ol <- .libPaths()
#' # called sequentially, .libPaths only add the last library
#' show( .libPaths('.') )
#' show( .libPaths(tempdir()) )
#' # restore
#' .libPaths(ol)
#' 
#' # .libPaths does not keep the standard user library
#' show( .libPaths() ) 
#' show( .libPaths('.') )
#' # restore
#' .libPaths(ol)
#' 
#' # with add_lib
#' show( add_lib('.') )
#' show( add_lib(tempdir()) )
#' show( add_lib('..', append=TRUE) )
#' 
#' # restore 
#' .libPaths(ol)
#' 
add_lib <- function(..., append=FALSE){
	
	p <- 
	if( append ) c(.libPaths(), ...)
	else c(..., .libPaths())
	.libPaths(p)
}


#' Package Check Utils
#' 
#' \code{isCRANcheck} \strong{tries} to identify if one is running CRAN-like checks.
#' 
#' Currently \code{isCRANcheck} returns \code{TRUE} if the check is run with 
#' either environment variable \code{_R_CHECK_TIMINGS_} (as set by flag \code{'--timings'})
#' or \code{_R_CHECK_CRAN_INCOMINGS_} (as set by flag \code{'--as-cran'}).
#' 
#' \strong{Warning:} the checks performed on CRAN check machines are on purpose not always 
#' run with such flags, so that users cannot effectively "trick" the checks.
#' As a result, there is no guarantee this function effectively identifies such checks.
#' If really needed for honest reasons, CRAN recommends users rely on custom dedicated environment 
#' variables to enable specific tests or examples.
#' 
#' @param ... each argument specifies a set of tests to do using an AND operator.
#' The final result tests if any of the test set is true.
#' Possible values are:
#' \describe{
#' \item{\code{'timing'}}{Check if the environment variable \code{_R_CHECK_TIMINGS_} is set, 
#' as with the flag \code{'--timing'} was set.}
#' \item{\code{'cran'}}{Check if the environment variable \code{_R_CHECK_CRAN_INCOMING_} is set, 
#' as with the flag \code{'--as-cran'} was set.}
#' }
#' 
#' @references Adapted from the function \code{CRAN}
#' in the \pkg{fda} package.
#' 
#' @export
isCRANcheck <- function(...){
  
  tests <- list(...)
  if( !length(tests) ){ #default tests
	  tests <- list('timing', 'cran')
  }
  test_sets <- c(timing="_R_CHECK_TIMINGS_", cran='_R_CHECK_CRAN_INCOMING_')
  tests <- sapply(tests, function(x){
			  # convert named tests
			  if( length(i <- which(x %in% names(test_sets))) ){
				  y <- test_sets[x[i]]
				  x <- x[-i]
				  x <- c(x, y)
			  }
			  # get environment variables
			  evar <- unlist(sapply(x, Sys.getenv))
			  all(nchar(as.character(evar)) > 0)
		  })
  
  any(tests)
}
#' \code{isCRAN_timing} tells if one is running CRAN check with flag \code{'--timing'}.
#' 
#' @export
#' @rdname isCRANcheck
isCRAN_timing <- function() isCRANcheck('timing')

#' \code{isCHECK} tries harder to test if running under \code{R CMD check}.
#' It will definitely identifies check runs for: 
#' \itemize{
#' \item unit tests that use the unified unit test framework defined by \pkg{pkgmaker} (see \code{\link{utest}});
#' \item examples that are run with option \code{R_CHECK_RUNNING_EXAMPLES_ = TRUE}, 
#' which is automatically set for man pages generated with a fork of \pkg{roxygen2} (see \emph{References}).
#' }
#' 
#' Currently, \code{isCHECK} checks both CRAN expected flags, the value of environment variable
#' \code{_R_CHECK_RUNNING_UTESTS_}, and the value of option \code{R_CHECK_RUNNING_EXAMPLES_}.
#' It will return \code{TRUE} if any of these environment variables is set to 
#' anything not equivalent to \code{FALSE}, or if the option is \code{TRUE}.
#' For example, the function \code{\link{utest}} sets it to the name of the package  
#' being checked (\code{_R_CHECK_RUNNING_UTESTS_=<pkgname>}), 
#' but unit tests run as part of unit tests vignettes are run with 
#' \code{_R_CHECK_RUNNING_UTESTS_=FALSE}, so that all tests are run and reported when 
#' generating them.
#' 
#' @references \url{https://github.com/renozao/roxygen2}
#' @rdname isCRANcheck
#' @export
#' 
#' @examples
#' 
#' isCHECK()
#' 
isCHECK <- function(){
	isCRANcheck() ||  # known CRAN check flags
            !isFALSE(utestCheckMode()) ||  # unit test-specific flag
            isTRUE(getOption('R_CHECK_RUNNING_EXAMPLES_')) # roxygen generated example flag
}

#' System Environment Variables
#' 
#' @param name variable name as a character string.
#' @param raw logical that indicates if one should return the raw value or
#' the convertion of any false value to \code{FALSE}.
#' 
#' @return the value of the environment variable as a character string or 
#' \code{NA} is the variable is not defined \strong{at all}.
#' 
#' @export
#' @examples
#' 
#' # undefined returns FALSE
#' Sys.getenv_value('TOTO')
#' # raw undefined returns NA
#' Sys.getenv_value('TOTO', raw = TRUE)
#' 
#' Sys.setenv(TOTO='bla')
#' Sys.getenv_value('TOTO')
#' 
#' # anything false-like returns FALSE
#' Sys.setenv(TOTO='false'); Sys.getenv_value('TOTO')
#' Sys.setenv(TOTO='0'); Sys.getenv_value('TOTO')
#' 
#' # cleanup
#' Sys.unsetenv('TOTO')
#' 
Sys.getenv_value <- function(name, raw = FALSE){
    val <- Sys.getenv(name, unset = NA, names = FALSE)
    if( raw ) return(val)
    # convert false values to FALSE if required
    if( is.na(val) || !nchar(val) || identical(tolower(val), 'false') || val == '0' ){
        val <- FALSE
    }
    val
}

checkMode_function <- function(varname){
    
    .varname <- varname
    function(value, raw = FALSE){
        if( missing(value) ) Sys.getenv_value(.varname, raw = raw)
        else{
            old <- Sys.getenv_value(.varname, raw = TRUE)
            if( is_NA(value) ) Sys.unsetenv(.varname) # unset
            else do.call(Sys.setenv, setNames(list(value), .varname)) # set value
            # return old value
            old	
        }	
    }
}


utestCheckMode <- checkMode_function('_R_CHECK_RUNNING_UTESTS_')
