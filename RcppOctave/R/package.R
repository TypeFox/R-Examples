#' @include utils.R
NULL

#' Seamless Interface to Octave -- and Matlab
#' 
#' The primary goal is to facilitate the port of Matlab/Octave scripts to R. 
#' The package enables to call any Octave functions from R and as well as browsing 
#' their documentation, passing variables between R and Octave, using R core RNGs 
#' in Octave, which ensure stochastic computations are also reproducible.
#'
#' \tabular{ll}{
#' Package: \tab RcppOctave\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2011-11-01\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @author
#' Renaud Gaujoux \email{renaud@@tx.technion.ac.il}
#'
#' Maintainer: Renaud Gaujoux \email{renaud@@tx.technion.ac.il}
#' @name RcppOctave
#' @rdname RcppOctave-package
#' @docType package
#' @keywords package
#' 
#' @cite Eaton2002
#' @bibliography ~/Documents/articles/library.bib
#' @examples
#' 
#' .CallOctave('svd', matrix(1:9, 3))
#' o_help('svd')
#' 
#' @import pkgmaker
#' @import Rcpp
#' @importFrom utils file_test
#' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
NULL

## #' @useDynLib RcppOctave

#inlineCxxPlugin <- function (...) 
#{
#	includes <- sprintf("%s\n#include <Rcpp.h>\n%s\n\n#ifndef BEGIN_RCPP\n#define BEGIN_RCPP\n#endif\n\n#ifndef END_RCPP\n#define END_RCPP\n#endif\n\nusing namespace Rcpp;\n", 
#			include.before, include.after)
#	list(env = list(PKG_LIBS = paste(libs, Rcpp:::RcppLdFlags())), 
#			includes = includes, LinkingTo = LinkingTo, body = function(x) {
#				sprintf("BEGIN_RCPP\n%s\nEND_RCPP", x)
#			}, Depends = Depends, Makevars = Makevars, Makevars.win = Makevars.win)
#}

# ##' Octave-RcppOctave Configuration
# ##' 
# ##' Configures Octave and load RcppOctave  
# ##' 
# ##' @param name Name of an RcppOctave path variable
# ##' @param ... extra names to be concatenated to the result with \code{\link{file.path}}.
# ##' Only used when \code{name} is not missing.
# ##' @param path path to Octave bin/ sub-directory
# ##' @return  a list (if \code{name is missing}) or a single character string.
# ##' 
# ##' @keywords internal
# ##' @export
# ##' @examples
# ##' 
# ##' OctaveConfig()
# ##' OctaveConfig('lib')
# ##' OctaveConfig('include')
# ##' OctaveConfig('modules')
# ##' 
#OctaveConfig <- local({
#	# config cache
#	.OctaveConfig <- NULL
#	function(name, ..., path=NULL){
#		    
#        do_init <- FALSE
#    	# return the whole config list if no name is provided
#    	if( is.null(.OctaveConfig) || !is.null(path) ){
#            
#            do_init <- TRUE
#            # save Octave bin path
#            if( !is.null(path) ) options(octave.path = path)
#                        
#            # create the config list at first call
#    		conf <- list(bin =  octave_config('BINDIR')
#                        , lib = octave_config(c('BINDIR', 'LIBDIR', 'OCTLIBDIR'))
#    					, include = octave_config(c('INCLUDEDIR', 'OCTINCLUDEDIR'))
#    			)
#    		
#    		# add a configuration variable for the module path
#    		conf$modules <- packagePath('modules')
#    		if( pkgmaker::isDevNamespace() ){ # fix module path due changes in devtools compilation step
#    			conf$modules <- file.path(tempdir(), packageName(), 'modules')
#    			# create module directory
#    			if( !file.exists(conf$modules) ){
#    				message("Faking devtools compilation directory '", conf$modules, "'")					
#    				dir.create(conf$modules, recursive=TRUE)
#    				src <- packagePath('src/modules')
#    				file.copy(file.path(src, c('PKG_ADD', list.files(src, pattern="*.oct$"))), conf$modules)
#    			}				
#    		} 
#    		
#            # save configuration
#    		.OctaveConfig <<- conf
#            
#    	}
#        
#        # setup/initialise library and modules if necessary
#        if( do_init ){
#            if( isTRUE(.OctaveInit()) )
#                .splash_message()
#        }
#        
#        if( missing(name) ){
#            if( !is.null(path) ) return( invisible(.OctaveConfig) )
#            else return( .OctaveConfig )
#        }
#    		
#    	settings <- .OctaveConfig[[name]]
#        # settings may contain more than one path => sapply
#    	unlist(sapply(settings, file.path, ...), use.names = FALSE)
#	}
#})

.debug <- function(...){
    if( any(!is.na(Sys.getenv(c('R_RCPPOCTAVE_VERBOSE', 'R_RCPPOCTAVE_DEBUG'), unset = NA))) ){ 
        message(...)
        TRUE
    }else FALSE
}

.load.lib <- function(libfile = NULL, pkgname, libname){
    
    lib_str <- 'RcppOctave library'
    if( !is.null(libfile) ) lib_str <- libfile
    
    if( .debug("  - Loading ", lib_str, " ... ", appendLF = FALSE) ){
        on.exit( .debug('ERROR') )
    }
    # load compiled library normally or in devmode
    res <- if( !is.null(libfile) ) dyn.load(libfile)
            else if( !isDevNamespace() ) library.dynam(pkgname, pkgname, libname)
        	else dyn.load(packagePath('src', paste0(pkgname, .Platform$dynlib.ext)))
    
    if( .debug('OK [', res[['path']], ']') ) on.exit()
    TRUE
}

# Load/Unload Octave Libraries
.OctaveLibs <- function(pkgname, libname){
    
    .load <- function(libfile = NULL){
        .load.lib(libfile, pkgname, libname)
    }
    
    .load_dep <- function( libdir = Octave.config[['libdir']] ){
        octlibs <- Octave.config[['libs']]
        dlibs <- file.first.path(libdir, paste0(octlibs, .Platform$dynlib.ext))
        if( length(notfound <- which(is.na(dlibs))) && is.Mac() ){
            # on Mac look for .dylib if default .so does not exist
            dlibs[notfound] <- file.first.path(libdir, paste0(octlibs[notfound], ".dylib"))
            notfound <- which(is.na(dlibs))
        }
		if( length(notfound) ){
			warning("RcppOctave - Could not find library ", paste0(octlibs[notfound], collapse = ','))
			dlibs <- dlibs[-notfound]
		}
        sapply(dlibs, .load)
    }
    
    dlls <- base::getLoadedDLLs()
	if ( pkgname %in%  names(dlls) ) return(TRUE)
    
    
    # custom installation
    if( isTRUE(Octave.config[['customed']]) ){
        .debug("* Explicit loading from custom location:")
        .load_dep()
        .load()
        return(TRUE)
    }
    
    # try directly loading the library
    .debug("* Try direct loading:")
    if( !is(try(.load(), silent = TRUE), 'try-error') ) return(TRUE)
    
    # check Octave configuration is reachable
    if( is.null(octave_bindir <- octave_config('BINDIR', mustWork = FALSE)) ){
        return(FALSE)
    }
    
    # setup path restoration for .onUnload
    .debug("* Expanding PATH with Octave lib directories:")
    on.exit( Sys.path$commit(), add = TRUE)
    Sys.path$append(octave_bindir)
    Sys.path$append(octave_config('OCTLIBDIR', mustWork = FALSE))
    Sys.path$append(octave_config('LIBDIR', exec = 'mkoctfile', mustWork = FALSE))
    .debug("  - ", Sys.path$get())
    
    # try reload
    .debug("* Re-try direct loading")
    if( !is(try(.load(), silent = TRUE), 'try-error') ) return(TRUE)
    
    # final try, pre-loading custom dependencies if necessary
    .debug("* Try excplicit loading:")
    .load_dep()
    .load()
    
    
    TRUE
}

.OctaveInit <- local({
    .ncall <- 0L
    .initialised <- FALSE
    function(libname = NULL, pkgname = NULL){
    
        # do not call recursively
        if( .ncall > 0L || .initialised ) return()
        .ncall <<- .ncall + 1L
        on.exit( .ncall <<- .ncall - 1L )
        #
        
        # force libname and pkgname when called from OctaveConfig
        if( is.null(pkgname) ) pkgname <- packageName() 
        if( is.null(libname) ) libname <- dirname(path.package(pkgname))
        #
        
    	# load required Octave libraries _before_ loading the package's library
    	if( !.OctaveLibs(pkgname, libname) ) return()
        
    	# start Octave session
    	octave_start(NULL)
        
    	# load Octave modules
    	octave_modules()
        
        # return TRUE
        .initialised <<- TRUE
    }
})

# dummy environment to trigger call to octave_end when quitting R
# via reg.finalizer (setup is done in .onLoad)
.octave_end_trigger <- environment()
.terminate_octave <- function(e){
    dlls <- base::getLoadedDLLs()
	if ( 'RcppOctave' %in%  names(dlls) ) octave_end()
    # unload dependencies
    sapply(Octave.config[['libs']], function(x){
                if( x %in%  names(dlls) ){
                    dyn.unload(dlls[[x]][['path']])
                }
            })
}

#' Platform Compatibility Check for RcppOctave
#' 
#' Checks if the current platform supports RcppOctave.
#' 
#' The following checks are performed:
#' \itemize{
#' \item platform must be either non-Windows or Windows-i386;  
#' \item The Octave version to be loaded is the same one as the one against
#' which RcppOctave libraries and modules were compiled.
#' }
#' 
#' @param details logical that indicates if the checks' details should be 
#' returned as well, in which case, the result is a list. 
#' 
#' @export
.isPlatformCompatible <- function(details = FALSE){
    
    res <- list(ok = TRUE, os.ok = TRUE, msg = NULL)
    .result <- function(res){
        if( details ) res
        else res[['ok']]
    }
    # check compatibility of OS
    res$os.ok <- .Platform$OS.type != 'windows' || .Platform$r_arch != 'x64'
    res$ok <- res$os.ok

    if( !res$ok && details ){
        pversion <- utils::packageVersion('RcppOctave')
        res$msg <- sprintf("RcppOctave [%s] - R platform %s is not supported", pversion, R.version$platform)
    }
    
    if( !res$ok ) return( .result(res) )
    
    # check compatibility of Octave compilation and loading version
    v <- octave_config('VERSION')
    v0 <- Octave.version[['version']]
    res$ok <- res$ok && utils::compareVersion(v0, v) == 0
    
    if( !res$ok && details ){
        res$msg <- c(res$msg, sprintf("RcppOctave was built against Octave %s, but loaded with Octave %s", v0, v))
    } 
    res$msg <- paste("NOTE:", res$msg, collapse = "\n")
    
    .result(res)
}

.onLoad <- function(libname, pkgname){

    # skip initialisation if not compatible platform
    compat <- .isPlatformCompatible(details = TRUE)
    if( !compat$ok ){
        .m <- message
        .m(compat$msg)
        # the library should load fine if the OS is not compatible 
        if( !compat$os.ok ) .load.lib(pkgname = pkgname, libname = libname)
        return()
    }
    
    # setup finalizer
    reg.finalizer(.octave_end_trigger, .terminate_octave, TRUE)
    
    # save initial PATH state to enable restoration in .onUnload
    Sys.path$init()
    
    # Initialise 
    .OctaveInit(libname, pkgname)
}

.onAttach <- function(libname, pkgname){
    
    .splash_message()
    
}

.splash_message <- function(){
    
    pversion <- utils::packageVersion('RcppOctave')
    
    compat <- .isPlatformCompatible(details = TRUE) 
    if( !compat$ok ){
        msg <- compat$msg
        packageStartupMessage(msg)
        if( interactive() ) warning(msg)
    
    }else if( is.null(octave_bindir <- octave_config('BINDIR', mustWork = FALSE)) ){
        packageStartupMessage("RcppOctave [", pversion, "] - Octave (not configured)"
                , "\nNOTE: Octave binaries were probably not found. See ?octave_config.")
    }else{
        # display info about config
        packageStartupMessage(sprintf("RcppOctave %s [Octave %s - path: %s]", pversion, o_version(), Octave.home('bin')))
    }
}

.onUnload <- function(libpath) {
	
    # skip cleanup if platform is not compatible
    if( !.isPlatformCompatible() ) return()
    
    # cleanup path on exit
    on.exit( Sys.path$revert("Reverting Octave changes to system PATH") )
    
    # terminate Octave session
    .terminate_octave(NULL)
    
	# unload compiled library normally or in devmode
	dlls <- base::getLoadedDLLs()
	pname <- 'RcppOctave'
	if ( pname %in%  names(dlls) ){
		if( !missing(libpath) )	library.dynam.unload(pname, libpath)
		else dyn.unload(dlls[[pname]][['path']])
	}
	
	# unload required Octave libraries
	#.OctaveLibs(FALSE)    
}

