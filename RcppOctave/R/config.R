# Project: RcppOctave
# 
# Octave configuration utils.
#
# Author: Renaud Gaujoux
# Created: Nov 25, 2013
###############################################################################

#' @include config-vars.R
NULL


#' Octave Configuration and Installation Information
#' 
#' The functions documented here enable retrieving information 
#' about the Octave installation used at installation or runtime
#' -- which should normally be the same.
#' 
#' @details \code{Octave.config} is a list that extends \code{Octave.version} with 
#' extra information about compilers and compilation flags. 
#' @rdname Octave.config
#' @export
#' 
#' @examples
#' Octave.config
Octave.config <- local({
    
    config <- setNames(.OCTAVE_CONFIG, tolower(names(.OCTAVE_CONFIG)))
    config$version.string = sprintf("Octave version %s (%s)", config$version, config$api_version)
    config$home <- dirname(config$bindir)
    config$libdir <- unique(c(config$libdir, Filter(nchar, gsub(" ", "", strsplit(config$lflags, "-L", fixed = TRUE)[[1L]]))))
    config$libs <- paste0('lib', Filter(nchar, gsub(" ", "", strsplit(config$libs, "-l", fixed = TRUE)[[1L]])))
    structure(config, class = 'simple.list')

})

#' @details \code{Octave.version} is list that contains version information as determined
#' by the configure script at installation time.
#'  
#' @rdname Octave.config
#' @family Octave.info
#' @export
#' 
#' @examples
#' Octave.version
Octave.version <- local({
    Octave.config[c('platform', grep('version', names(Octave.config), fixed = TRUE, value = TRUE))]
})

#' Octave Home Directory
#' 
#' Returns the path to Octave home directory, i.e. the directory that contains
#' the \code{bin/} sub-directory where Octave binaries can be found, 
#' e.g., typically \code{/usr/} on Linux machines. 
#' 
#' The path to Octave home directory is determined in the following 
#' order:
#' 
#' \itemize{
#' \item value of global option \code{'octave.home'};
#' \item value of the environment variable \code{'OCTAVE_HOME'};
#' \item path used during configuration/installation of \pkg{RcppOctave}.
#' \item path returned by \code{octave-config}, which is itself looked up 
#' in the system PATH.
#' }
#' 
#' If set, the global option or environment variable should contain 
#' the \strong{absolute} path to Octave root directory.
#' 
#' @param ... character strings that are appended to Octave path
#' via \code{\link{file.path}}.
#' @param configure logical that indicates if one should directly  
#' return the path that was used when configuring (i.e. installing) 
#' \pkg{RcppOctave}
#' @param use.system logical that indicates if one should try
#' using \code{octave-config} to retrieve Octave home directory.
#' This would be done as a last resort, if the path could not be 
#' determined in any other ways (see \emph{Details}).
#' 
#' @return a character string, or \code{NULL} if the path was not found.
#' 
#' @family Octave.info
#' @export 
#' @examples 
#' Octave.home()
#' 
Octave.home <- function(..., configure = Octave.config[['customed']], use.system = TRUE){
    .octave_home_confvar <- sprintf('@%s@', 'OCTAVE_BINDIR')

    # default is to use the path used at installation
    .config.path <- Octave.config[['home']]
    if( identical(.config.path, .octave_home_confvar) ) .config.path <- ''
    if( configure ) return(file.path(.config.path, ...))
    
    # check global R option
    if( is.null(path <- getOption('octave.home')) ){
        # check environment variable OCTAVE_HOME
        if( !nzchar(path <- Sys.getenv('OCTAVE_HOME')) ){
            # check existence of path resolved at configure time
            path <- if( Octave.config[['customed']] ) .config.path else ''
            if( (nzchar(path) || !file_test('-d', path)) && use.system ){
                # last resort: retrieve path from octave-config
                path <- dirname(octave_config('BINDIR', mustWork = FALSE, warn = FALSE, bindir = NA))
            }
        }
    }
    
    if( length(path) && nzchar(path) ) file.path(path, ...) 
}


modules.path <- function(){
    
    modpath <- packagePath('modules')
    # fix module path due changes in devtools compilation step
    if( isDevNamespace() ) modpath <- file.path(tempdir(), packageName(), 'modules')
    # append architecture sub-directory if necessary
    if( nzchar(.Platform$r_arch) ) modpath <- file.path(modpath, .Platform$r_arch)
    modpath
}

#' Octave Session Details
#' 
#' \code{Octave.info} is a function that retrieves information about 
#' the version of Octave that is used by the current session 
#' of \pkg{RcppOctave}. 
#' 
#' @param name name of the variable to retrieve
#' 
#' @rdname Octave.config
#' @export
#' @examples 
#' Octave.info()
#' Octave.info('modules')
#' 
Octave.info <- function(name){
    
    # special handling of module path
    modpath <- modules.path()
    if( isDevNamespace() ){
        # create module directory
        if( !file.exists(modpath) ){
            message("Faking devtools compilation directory '", modpath, "'")					
            dir.create(modpath, recursive=TRUE)
            src <- packagePath('src/modules')
            file.copy(file.path(src, c('PKG_ADD', list.files(src, pattern="*.oct$"))), modpath)
        }				
    }

    res <- list(
                version = o_version(),
                home = Octave.home(),
                home.config = Octave.home(configure = TRUE),
                modules = modpath
                )
    
    # add some Octave configuration variables
    live_info <- o_config_info(c(cc = 'CC', cc.version = 'CC_VERSION', f77 = 'F77'))
    res <- structure(c(res, as.list(live_info)), class = 'simple.list')
    if( !missing(name) ){
        res <- res[[name]]
    }
    res
}


