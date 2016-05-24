#' Determine SAGA GIS default paths
#'
#' Internal functions that determine OS-specific default paths in which SAGA GIS binaries and modules might be located.
#' @name rsaga.default.path
#' @param sysname character: name of the operating system, determined by default by \code{\link[base]{Sys.info}}: e.g., \code{"Windows"}, \code{"Linux"}, \code{"Darwin"} (for Mac OSX), or \code{"FreeBSD"}
#' @details These functions are used internally by \code{\link{rsaga.env}}. On Windows systems, the guess made by \code{rsaga.default.path} is \code{"C:/Progra~1/SAGA-GIS"}, and for the modules, it is the \code{"modules"} subfolder of the path with the binaries.
#' 
#' On non-Windows systems, \code{rsaga.default.path} submits a \code{which saga_cmd} call to the operating system to find the binaries, usually in \code{/usr/local/bin} or in \code{/usr/bin}. To find the modules, \code{rsaga.default.modules.path} first checks if a \code{SAGA_MLB} environment variable exists. If not, it will replace the \code{/bin} part (if present) with \code{/lib/saga} or otherwise it just guesses that it's \code{/usr/local/lib/saga}.
#' @seealso \code{\link{rsaga.env}}
#' @export
rsaga.default.path = function(sysname = Sys.info()["sysname"]) {
    if (sysname == "Windows") {
        res = "C:/Progra~1/SAGA-GIS"
    } else { ### tested with: ((sysname == "Linux") | (sysname == "Darwin") | (sysname == "FreeBSD"))
        res = sub('/saga_cmd','',system2('which',args='saga_cmd',stdout=TRUE))
    }
    ###} else stop("operating system '", sysname, "' not recognized or not supported\n",
    ###        "try specifying SAGA GIS path and module path manually")
    return(res)
}

#' @rdname rsaga.default.path
#' @name rsaga.default.modules.path
#' @param saga.path character: path with SAGA GIS binaries, as determined (e.g.) by \code{rsaga.default.path}
#' @export
rsaga.default.modules.path = function(sysname = Sys.info()["sysname"], 
    saga.path = rsaga.default.path(sysname))
{
    if (sysname == "Windows") {
        modules = file.path(saga.path,"modules")
    } else { ### tested with: ((sysname == "Linux") | (sysname == "Darwin") | (sysname == "FreeBSD"))
        modules = Sys.getenv("SAGA_MLB")[[1]]
        if (modules == "") {
            # have a backup path in case SAGA_MLB is not set/empty
            if (substr(saga.path, nchar(saga.path)-3, nchar(saga.path)) == paste(.Platform$file.sep, "bin", sep="")) {
                modules = file.path( substr(saga.path, 1, nchar(saga.path)-4), "lib", "saga" )
            } else modules = "usr/local/lib/saga"
        }
    }
    ###} else stop("operating system '", sysname, "' not recognized or not supported\n",
    ###        "try specifying SAGA GIS path and module path manually")
    return(modules)
}



#' Set up the RSAGA Geoprocessing Environment
#'
#' \code{rsaga.env} creates a list with system-dependent information on SAGA path, module path and data  (working) directory. This kind of a list is required by most RSAGA geoprocessing functions and is referred to as the 'RSAGA geoprocessing environment.'
#' @name rsaga.env
#' @param workspace path of the working directory for SAGA; defaults to the current directory (\code{"."}).
#' @param cmd name of the SAGA command line program; defaults to \code{saga_cmd.exe}, its name under Windows
#' @param path path in which to find \code{cmd}; \code{rsaga.env} is usually able to find SAGA on your system if it is installed; see Details.
#' @param modules path in which to find SAGA libraries; see Details
#' @param version optional character string: SAGA GIS (API) version, e.g. \code{"2.0.8"}; if missing, a call to \code{\link{rsaga.get.version}} is used to determine version number of SAGA API
#' @param cores optional numeric argument, or \code{NA}: number of cores used by SAGA GIS; supported only by SAGA GIS 2.1.0 (and higher), ignored otherwise (with a warning). Multicore-enabled SAGA GIS modules such as the one used by \code{\link{rsaga.pisr}} seem to run in multicore mode by default when this argument is not specified, therefore \code{cores} should only be specified to use a smaller number of cores than available on a machine.
#' @param parallel optional logical argument (default: \code{FALSE}): if \code{TRUE}, run RSAGA functions that are capable of parallel processing in parallel mode; note that this is completely independent of the behaviour of SAGA GIS (which can be controlled using the \code{cores} argument); currently only some RSAGA functions support parallel processing (e.g., \code{\link{pick.from.ascii.grid}} or \code{\link{rsaga.get.modules}}). \code{parallel=TRUE} requires that a parallel backend such as \pkg{doSNOW} or \pkg{doMC} is available and has been started prior to calling any parallelized RSAGA function, otherwise warnings may be generated
#' @param check.libpath if \code{TRUE} (default), first look for SAGA GIS in the folder where the RSAGA package is installed
#' @param check.SAGA if \code{TRUE} (default), next check the path given by the environment variable \code{SAGA}, if it exists
#' @param check.PATH if \code{TRUE} (default on Windows), next look for SAGA GIS in all the paths in the \code{PATH} environment variable; defaults to \code{FALSE} on non-Windows OS
#' @param check.os.default if \code{TRUE}, look for SAGA GIS in the folder specified by \code{os.default.path}.
#' @param os.default.path on Windows, \code{C:/Progra~1/SAGA-GIS}; on unix, an attempt is made to locate \code{saga_cmd}
#' @param lib.prefix character string: a possible (platform-dependent) prefix for SAGA GIS library names; if missing (recommended), a call to \code{\link{rsaga.lib.prefix}} tries to determine the correct prefix, e.g. \code{""} on Windows systems and \code{"lib"} on non-Windows systems with SAGA GIS pre-2.1.0. Try specifying \code{""} or \code{"lib"} manually if this causes problems, and contact the package maintainer if the detection mechanism fails on your system (indicate your \code{Sys.info()["sysname"]} and your SAGA GIS version)
#' @details IMPORTANT: Unlike R functions such as \code{\link{options}},  which changes and saves settings somewhere in a global variable, \code{\link{rsaga.env}} does not actually 'save' any settings, it simply creates a list that can (and has to) be passed to other \code{rsaga.*} functions. See example below.
#' 
#' I strongly recommend to install SAGA GIS in \code{"C:/Program Files/SAGA-GIS"}  in the case of English-language Windows platforms (the equivalent non-English installation folder in the case of non-English Windows versions seems to work as well). If this is the only SAGA GIS copy on the computer and you do \emph{not} define a Windows environment variable \code{SAGA}, then RSAGA should normally be able to find your SAGA GIS installation in this folder.
#' 
#' \code{rsaga.env} tries to collect infromation on the (R)SAGA environment. If \code{path} is missing, \code{rsaga.env} first looks for an environment variable \code{SAGA}; if this is undefined, it checks the current working directory, then the paths given in the PATH environment variable, and finally the function's guess is \code{"C:/Progra~1/SAGA-GIS"} (or \code{"/usr/local/bin"} on non-Windows systems).
#' 
#' The default \code{modules} folder on Windows systems is the  \code{modules} subfolder of the SAGA binaries' folder. The \code{SAGA_MLB} environment variable is \emph{not} checked by \code{rsaga.env}.
#'
#' On Unix (and Mac OS X) systems, the default \code{modules} folder is as specified in the \code{SAGA_MLB} environment variable. If this is empty / not set, then the following backup path is used. If \code{path} ends with "/bin", then "/bin" is changed to "/lib/saga" and taken as the \code{modules} path; otherwise, \code{/usr/local/lib/saga} is used.
#' @return A list with components \code{workspace}, \code{cmd}, \code{path}, \code{modules}, \code{version}, \code{cores} and \code{parallel} with values as passed to \code{rsaga.env} or default values as described in the Details section.
#' @author Alexander Brenning
#' @note Note that the default \code{workspace} is \code{"."}, not \code{getwd()}; i.e. the default SAGA workspace folder is not fixed, it changes each time you change the R working directory using \code{setwd}.
#' @seealso \code{\link{rsaga.get.version}}
#' @examples
#' \dontrun{
#' # Check the default RSAGA environment on your computer:
#' myenv <- rsaga.env()
#' myenv
#' # SAGA data in C:/sagadata, binaries in C:/SAGA-GIS, modules in C:/SAGA-GIS/modules:
#' myenv <- rsaga.env(workspace="C:/sagadata", path="C:/SAGA-GIS")
#' # Unix: SAGA in /usr/bin (instead of the default /usr/local/bin),
#' # and modules in /use/lib/saga:
#' # myenv <- rsaga.env(path="/usr/bin")
#' # Use the 'myenv' environment for SAGA geoprocessing:
#' rsaga.hillshade("dem","hillshade",env=myenv)
#' # ...creates (or overwrites) grid "C:/sagadata/hillshade.sgrd"
#' # derived from digital elevation model "C:/sagadata/dem.sgrd"
#' 
#' # Same calculation with different SAGA version:
#' # (I keep several versions in SAGA-GIS_2.0.x folders:)
#' myenv05 = rsaga.env(path = "C:/Progra~1/SAGA-GIS_2.0.5")
#' rsaga.hillshade("dem","hillshade205",env=myenv05)
#' }
#' @keywords spatial interface
#' @export
rsaga.env = function( workspace=".", 
    cmd = ifelse(Sys.info()["sysname"]=="Windows", "saga_cmd.exe", "saga_cmd"), 
    path, modules, version, cores, parallel = FALSE,
    check.libpath = TRUE, 
    check.SAGA = TRUE, 
    check.PATH = Sys.info()["sysname"]=="Windows",
    check.os.default = TRUE,
    os.default.path = rsaga.default.path(),
    lib.prefix )
    #os.default.path = ifelse(.Platform$OS.type=="windows",
    #    "C:/Progra~1/SAGA-GIS",
    #    sub('/saga_cmd','',system2('which',args='saga_cmd',stdout=TRUE))) )
{
    rsaga.get.possible.SAGA.paths = function( check.libpath, check.SAGA,
            check.PATH, check.os.default, os.default.path )
    {
        path = c()
        if (check.libpath) {
            path.pckg = file.path( .libPaths(), "RSAGA" )
            try( path.pckg <- path.package("RSAGA"), silent = TRUE )
            path = c( path, file.path( path.pckg, "saga_vc" ) )
            path = c( path, file.path( path.pckg, "SAGA-GIS" ) )
        }
        if (check.SAGA)
            path = c( path, Sys.getenv("SAGA"), getwd() )

        if (check.PATH)
        {
            os.path = Sys.getenv("PATH")
    
            # This won't work:
            #os.paths = strsplit(os.path,";")[[1]]
            # on Windows, filenames may contain a ";"!
    
            # Split the PATH environment variable at ";",
            # except when the ";" occurs inside a path in quotation marks:
            i = 1
            while (i <= nchar(os.path)) {
                next.path = ""
                while ((i <= nchar(os.path)) & (substr(os.path,i,i) != ";")) {
                    # repeat until end of PATH or ";":
                    if (substr(os.path,i,i) == '"') { ## potential problems on unix??
                        # if a quote begins, ignore any ";"s:
                        i = i + 1 # skip the quotation mark
                        while ((i <= nchar(os.path)) & (substr(os.path,i,i) != '"')) {
                            # inside the quotation, pass any normal characters
                            # and possible ";'s to next.path:
                            next.path = paste(next.path, substr(os.path,i,i), sep="")
                            i = i + 1
                        } # until end of PATH or end of quotation
                        if (substr(os.path,i,i) != '"') { ## potential problems on unix??
                            # unlikely case that quotation ends at the end of PATH,
                            # with the closing quotation mark missing:
                            next.path = paste(next.path, substr(os.path,i,i), sep="")
                        }
                        # note that the '"' is not passed to next.path
                    } else {
                        # pass any 'normal' character from PATH to next.path:
                        next.path = paste(next.path, substr(os.path,i,i), sep="")
                    }
                    i = i + 1
                } # until end of PATH or next ";"
                path = c(path, next.path)
                i = i + 1
            } # until end of PATH
        }
        
        if (check.os.default) {
            if (!is.null(os.default.path))
                path = c(path, os.default.path)
        }

        return(path)
    }
    
    
    # Default RSAGA workspace is the current working directory of R:
    if (workspace == "") workspace = "."
    
    # cmd has inheritted a 'sysname' name from the result of Sys.info():
    cmd = unname(cmd)

    # If path specified by user, check if valid:
    if (missing(path)) {
        path = NULL
    } else if (!file.exists(file.path(path,cmd))) {
        warning("SAGA command line program ", cmd, " not found in the specified path ", 
            path, ".", "\nTrying to find it somewhere else.")
        path = NULL
    }
    # No path information available, start search...
    if (is.null(path)) {
        # Ordered list of candidate paths:
        path = rsaga.get.possible.SAGA.paths( check.libpath, check.SAGA,
                                    check.PATH, check.os.default,
                                    os.default.path = os.default.path )
        if (length(path)==0) {
            warning("don't know where to look for SAGA command line program, no paths specified",
                "\nTrying to find it in the current working directory...")
            path = getwd()
        }
        # Determine first candidate folder that contains saga_cmd executable:
        the.path = NULL
        for (pa in path) {
            if (file.exists(file.path(pa,cmd))) {
                the.path = pa
                break
            }
        }
        # None found:
        if (is.null(the.path)) {
            warning("SAGA command line program '", cmd, 
                    "' not found in any of the following paths\n", paste(path,collapse="\n"))
            return(NULL)
        }
        path = the.path
    }    

    # Set up module path:
    if (missing(modules)) {
        modules = rsaga.default.modules.path(saga.path = path)
    } else {
        # Empty character string interpreted as ".",
        # i.e. current working directory:
        if (modules == "") modules = getwd()
    }
    
    # Check if folder with the specified names exist:
    if (!file.exists(workspace))
        warning("Invalid workspace path ", workspace)
    if (!file.exists(path))
        warning("Invalid SAGA path ", path)
    if (!file.exists(modules))
        warning("Invalid SAGA modules path ", modules)

    ## Check if saga_cmd[.exe] exists in specified folder:
    ## (no need to do this - this was done above)
    #if (!file.exists(file.path(path,cmd)))
    #    warning("SAGA GIS command line program not found.\nFile name: ", file.path(path,cmd))

    # (optional) number of cores:
    if (missing(cores)) {
        cores = NA
    } else stopifnot(cores >= 0)

    # Set up RSAGA geoprocessing environment:
    env = list(
        workspace = workspace,
        cmd = cmd,
        path = path,
        modules = modules,
        version = NA,
        cores = cores,
        parallel = parallel )
        
    # Determine SAGA API version, if not specified by caller:
    if (missing(version))
        version = rsaga.get.version(env = env)
    env$version = version
    
    if (!is.na(env$cores) & (is.na(env$version) | (substr(env$version,1,4) == "2.0."))) {
        warning("'cores' argument not supported by SAGA GIS versions <2.1.0; ignoring 'cores' argument\n",
            "Use SAGA GIS 2.1.0+ for multicore geoprocessing.")
        env$cores = NA
    }

    if (missing(lib.prefix))
        lib.prefix = rsaga.lib.prefix(env = env)
    env$lib.prefix = lib.prefix
    
    return( env )
}


#' Determine prefix for SAGA GIS library names
#'
#' Internal function that determines the possible prefix for SAGA GIS library names - relevant for non-Windows SAGA GIS pre-2.1.0.
#' 
#' @name rsaga.lib.prefix
#' @param env list, setting up a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}.
#' @details Some non-Windows versions of \code{saga_cmd} require library names with a \code{"lib"} prefix, e.g. \code{libio_grid} instead of \code{io_grid}. This function, which is called by \code{\link{rsaga.env}} tries to guess this behaviour based on the operating system and SAGA GIS version.
#' @return A character string, either \code{""} or \code{"lib"}.
#' @seealso \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' env = rsaga.env()
#' # obtained by a call to rsaga.lib.prefix:
#' env$lib.prefix
#'
#' # more explicitly:
#' rsaga.lib.prefix(env=env)
#' }
#' @keywords spatial interface
#' @export
rsaga.lib.prefix = function(env) {
    lib.prefix = "lib"
    if ((Sys.info()["sysname"] == "Windows")) {
        lib.prefix = ""
    } else if ((Sys.info()["sysname"] == "Darwin")) {
        lib.prefix = ""
    } else if (!is.na(env$version)) {
        if (substr(env$version,1,4) == "2.1." | substr(env$version,1,4) == "2.2.")
            lib.prefix = ""
    }
    return(lib.prefix)
}




#' Determine SAGA GIS version
#'
#' Determine SAGA GIS version.
#' 
#' @name rsaga.get.version
#' @param env list, setting up a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}. Note that \code{version=NA} ensures that \code{\link{rsaga.env}} won't call \code{rsaga.get.version} itself.
#' @param ... additional arguments to \code{\link{rsaga.geoprocessor}}
#' @details The function first attempts to determine the SAGA version directly through a system call \code{saga_cmd --version}, which is supported by SAGA GIS 2.0.8+. If this fails, \code{saga_cmd -h} is called, and it is attempted to extract the version number of the SAGA API from the output generated, which works for 2.0.4 - 2.0.7.
#' @return A character string defining the SAGA GIS (API) version. E.g., \code{"2.0.8"}.
#' @seealso \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' myenv <- rsaga.env()
#' myenv$version
#' # rsaga.env actually calls rsaga.get.version:
#' rsaga.get.version()
#' 
#' # I keep several versions of SAGA GIS in SAGA-GIS_2.0.x folders:
#' myenv05 = rsaga.env(path = "C:/Progra~1/SAGA-GIS_2.0.5", version = NA)
#' # Check if it's really version 2.0.5 as suggested by the folder name:
#' rsaga.get.version(env = myenv05)
#' }
#' @keywords spatial interface
#' @export
rsaga.get.version = function(env = rsaga.env(version=NA), ...) 
{
    version = NA
    
    reduce.version = function(x) {
        x = gsub(" ", "", x)
        # reduce e.g. 2.1.0(32bit) to 2.1.0; added 2013-05-28
        x = strsplit(x, '(', fixed=TRUE)[[1]][1]
        return(x)
    }

    # Added 27-Dec-2011:
    # saga_cmd --version (only works in SAGA GIS 2.0.8+)
    out = rsaga.geoprocessor(lib = NULL, prefix = "--version", show.output.on.console = FALSE, 
        flags = NULL, warn = -1, env = env, ...)
    if (all(out != "error: module library not found [--version]")) {
        if (length(out >= 1)) {
            if (any(sel <- (substr(out,1,9) == "SAGA API:"))) {
                # This first option was mentioned by Volker Wichmann on [saga-gis-developer]
                # although SAGA GIS 2.0.8 for Windows uses a different output format:
                out = gsub("SAGA API:", "", out[sel][1], fixed = TRUE)
                out = reduce.version(out)
                return(out)
            } else if (any(sel <- (substr(out,1,13) == "SAGA Version:"))) {
                # Output format used by SAGA GIS 2.0.8 for Windows:
                # SAGA Version: 2.0.8
                out = gsub("SAGA Version:", "", out[sel][1], fixed = TRUE)
                out = reduce.version(out)
                return(out)
            }
        }
    }
    # End added code

    # Older SAGA GIS versions:
    # ------------------------
    
    #### check if this function works on unix????
    # Retrieve basic help page of saga_cmd:
    out = rsaga.geoprocessor(lib = NULL, prefix = "-h", show.output.on.console = FALSE, 
        flags = NULL, warn = -1, env = env, ...)
        
    # Process the help page line by line in order to find lines starting
    # with "SAGA API " (or "SAGA CMD ", if no "SAGA API " line available)
    for (i in 1:length(out)) {
        if (substr(out[i],1,9) == "SAGA API ") {
            if (as.numeric(substr(out[i],10,10)) > 0) {
                version = gsub(" ", "", substr(out[i],10,nchar(out[i])), fixed = TRUE)
                break
            }
        } else if (substr(out[i],1,9) == "SAGA CMD ") {
            if ( is.na(version) & any( as.character(0:9) == substr(out[i],10,10) ) ) {
                version = gsub(" ", "", substr(out[i],10,nchar(out[i])), fixed = TRUE)
                # no 'break' here because we're still hoping to find info on SAGA API version
                # SAGA 2.0.4 only shows SAGA CMD version = 2.0.4, however *some* versions
                # distinguish between SAGA API and SAGA CMD versions.
            }
        }
    }
    return(version)
}



#' Find SAGA libraries and modules
#'
#' These functions list the SAGA libraries (\code{rsaga.get.libraries}) and modules (\code{rsaga.get.lib.modules}, \code{rsaga.get.modules}) available in a SAGA installation, and allow to perform a full-text search among these functions.
#' @name rsaga.get.modules
#' @param text character string to be searched for in the names of available libraries and/or modules
#' @param search.libs logical (default \code{TRUE}); see \code{search.modules}
#' @param search.modules logical (default \code{TRUE}): should \code{text} be searched for in library and/or module names?
#' @param ignore.case logical (default \code{FALSE}): should the text search in library/module names be case sensitive?
#' @param lib character string with the name of the library in which to look for modules
#' @param libs character vector with the names of libraries in which to look for modules; if missing, all libraries will be processed
#' @param module module name or numeric code
#' @param modules optional list: result of \code{rsaga.get.modules}; if missing, a list of available modules will be retrieved using that function
#' @param env a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}
#' @param path path of SAGA library files (\code{modules} subfolder in the SAGA installation folder); defaults to the path determined by \code{\link{rsaga.env}}.
#' @param dll file extension of dynamic link libraries
#' @param interactive logical (default \code{FALSE}): should modules be returned that can only be executed in interactive mode (i.e. using SAGA GUI)?
#' @param parallel logical (defaults to \code{env$parallel}): if \code{TRUE}, run in parallel mode; requires a parallel backend such as \pkg{doSNOW} or \pkg{doMC}
#' @param ... currently only \code{interactive} to be passed on to \code{rsaga.get.lib.modules}
#' @return \code{rsaga.get.libraries} returns a character vector with the names of all SAGA libraries available in the folder \code{env$modules}.
#'
#' \code{rsaga.get.lib.modules} returns a \code{data.frame} with:
#' \itemize{
#' \item{name} {the names of all modules in library \code{lib},}
#' \item{code} {their numeric identifiers,}
#' \item{interactive} {and a logical variable indicating whether a module can only be executed in interactive (SAGA GUI) mode.}
#' }
#'
#' \code{rsaga.get.modules} returns a list with, for each SAGA library in \code{libs}, a \code{data.frame} with module information as given by \code{rsaga.get.lib.modules}. If \code{libs} is missing, all modules in all libraries will be retrieved.
#' 
#' @note For information on the usage of SAGA command line modules, see \code{\link{rsaga.get.usage}}, or \code{\link{rsaga.html.help}} (in SAGA GIS 2.1.0+), or the RSAGA interface function, if available.
#' @seealso \code{\link{rsaga.get.usage}}, \code{\link{rsaga.html.help}}, \code{\link{rsaga.geoprocessor}}, \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' # make sure that 'rsaga.env' can find 'saga_cmd.exe'
#' # before running this:
#' rsaga.get.libraries()
#' # list all modules in my favorite libraries:
#' rsaga.get.modules(c("io_grid", "grid_tools", "ta_preprocessor",
#'     "ta_morphometry", "ta_lighting", "ta_hydrology"))
#' # list *all* modules (quite a few!):
#' # rsaga.get.modules(interactive=TRUE)
#' 
#' # find modules that remove sink from DEMs:
#' rsaga.search.modules("sink")
#' # find modules that close gaps (no-data areas) in grids:
#' rsaga.search.modules("gap")
#' }
#' @keywords spatial interface
#' @export
rsaga.get.modules = function(libs, env = rsaga.env(), 
    interactive = FALSE, parallel = env$parallel)
{
    if (missing(libs)) libs = rsaga.get.libraries(env$modules)
    op = options(warn = -1) # llply would generate two warnings
    on.exit(options(op))    
    res = llply(libs, .fun=rsaga.get.lib.modules, env = env, 
        interactive = interactive, .parallel = parallel)
    names(res) = libs
    return(res)
}

#' @rdname rsaga.get.modules
#' @name rsaga.get.libraries
#' @export
rsaga.get.libraries = function(path = rsaga.env()$modules, dll)
{
    if (missing(dll)) {
        dll = .Platform$dynlib.ext
        if (Sys.info()["sysname"] == "Darwin") dll = ".dylib"
    }
    dllnames = dir(path,paste("^.*\\",dll,"$",sep=""))
    if (Sys.info()["sysname"] != "Windows") ### %in% c("Linux","Darwin","FreeBSD"))
        if (all(substr(dllnames,1,3) == "lib"))
            dllnames = substr(dllnames, 4, nchar(dllnames)) # remove the "lib"
    return( gsub(dll,"",dllnames,fixed=TRUE ) )
}


#' @rdname rsaga.get.modules
#' @name rsaga.get.lib.modules
#' @export
rsaga.get.lib.modules = function(lib, env=rsaga.env(), interactive=FALSE)
{
    res = NULL

    # changed by Rainer Hurling, 2013-07-23:
    ###if ( lib == "opencv" & (is.na(env$version) | (env$version == "2.0.4" | env$version == "2.0.5" | env$version == "2.0.6")) ) {
    if ( lib == "opencv" & env$version %in% c(NA,"2.0.4","2.0.5","2.0.6") ) {
        warning("skipping library 'opencv' because it produces an error\n",
            "  when requesting its module listing in SAGA version 2.0.4 - 2.0.6)")
        # return an empty data.frame of the same format as in the successful situation:
        return( data.frame( code = numeric(), name = character(), interactive = logical() ) )
    }

    rawres = rsaga.geoprocessor(lib, module=NULL, env=env,
        intern=TRUE, show.output.on.console=FALSE, flags=NULL, invisible=TRUE,
        reduce.intern=FALSE, check.module.exists=FALSE, warn = -1)

    wh = which( gsub(" ","",tolower(rawres)) %in% c("availablemodules:","executablemodules:","modules:", "tools:") )

    if (length(wh) > 0) {
        rawres = rawres[ (wh[length(wh)]+1) : length(rawres) ]
        rawres = rawres[ rawres != "" ]
        rawres = rawres[ rawres != "type -h or --help for further information" ]
        # inserted tolower() for SAGA 2.1.0 RC1:
        rawres = rawres[ tolower(rawres) != "error: module" ]
        rawres = rawres[ tolower(rawres) != "error: tool" ]
        rawres = rawres[ tolower(rawres) != "error: select a tool" ]
    }
    if (length(wh) > 0) {
        rawres = strsplit(rawres,"\t- ")
        mcodes = c()
        mnames = c()
        minteracs = c()
        for (descr in rawres) {
            mygrep = c( grep("[",descr[1],fixed=TRUE), grep("]",descr[1],fixed=TRUE),
                grep("[interactive]",descr[2],fixed=TRUE) )
            minterac = (length(mygrep) > 0)
            # skip interactive modules if only interactive ones are allowed:
            if (!minterac | interactive) {
                mcode = gsub("[","",gsub("]","",gsub(" ","",descr[1]),fixed=TRUE),fixed=TRUE)
                mname = gsub("[interactive] ","",descr[2],fixed=TRUE)
                mcodes = c(mcodes, as.numeric(mcode))
                mnames = c(mnames, mname)
                minteracs = c(minteracs, minterac)
            }
        }
        #if (length(mcodes) > 0)
        res = data.frame(code=mcodes, name=mnames, interactive=minteracs)
    }
    return(res)
}

#' @rdname rsaga.get.modules
#' @name rsaga.module.exists
#' @export
rsaga.module.exists = function(libs, module, env = rsaga.env(), ...) {
    if (missing(libs)) libs = rsaga.get.libraries(env$modules)
    wh = "name"
    if (is.numeric(module)) wh = "code"
    for (i in 1:length(libs)) {
        modules = rsaga.get.lib.modules(libs[i], env = env, ...)
        if (!is.null(modules))
            if (any(modules[,wh] == module))
                return(TRUE)
    }
    return(FALSE)
}


#' @rdname rsaga.get.modules
#' @name rsaga.search.modules
#' @export
rsaga.search.modules = function(text, modules, search.libs=TRUE, search.modules=TRUE,
    env=rsaga.env(), ignore.case=TRUE, ...)
{
    pattern = paste("^.*",text,sep="")
    lib = NULL
    mod = NULL
    if (search.libs) {
        lib.nm = rsaga.get.libraries(path=env$modules)
        wh.lib = grep(pattern,lib.nm,ignore.case=ignore.case)
        lib = lib.nm[wh.lib]
    }
    if (search.modules) {
        if (missing(modules))  
            modules = rsaga.get.modules(env=env,...)
        mod.nm = unlist(sapply(modules,function(x) if (is.atomic(x)) NULL else as.character(x$name)),use.names=FALSE)
        mod.libs = sapply(modules,function(x) if (is.atomic(x)) 0 else nrow(x))
        mod.libs = rep(names(mod.libs),mod.libs)
        wh.mod = grep(pattern,mod.nm,ignore.case=ignore.case)
        mod = data.frame( lib=mod.libs[wh.mod], module=mod.nm[wh.mod] )
    }
    return( list( lib = lib, modules = mod ) )
}




#' Usage of SAGA command line modules
#' 
#' \code{rsaga.get.usage} provides information on the usage of and arguments required by SAGA command line modules.
#'
#' @name rsaga.get.usage
#' @param lib name of the SAGA library
#' @param module name or numeric identifier of SAGA module in library \code{lib}
#' @param env a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}
#' @param show logical (default: \code{TRUE}); display usage in the R console?
#'
#' @details This function is intended to provide information required to use the 
#' \code{\link{rsaga.geoprocessor}} and for writing your own high-level interface 
#' function for SAGA modules. R--SAGA interfaces already exist for some SAGA modules, 
#' e.g. \code{\link{rsaga.hillshade}}, \code{\link{rsaga.local.morphometry}}, but there 
#' are many more.
#' @return The character vector with usage information is invisibly returned.
#' @seealso \code{\link{rsaga.html.help}}, \code{\link{rsaga.geoprocessor}}, \code{\link{rsaga.env}}, \code{\link{rsaga.get.modules}}
#' @examples
#' \dontrun{
#' rsaga.get.usage("io_grid",1)
#' rsaga.get.usage("ta_preprocessor",2)
#' rsaga.get.usage("ta_morphometry",0)
#' # in SAGA GIS 2.1.0+, compare:
#' rsaga.html.help("io_grid",1)
#' # etc.
#' }
#' @keywords spatial interface
#' @export
rsaga.get.usage = function(lib, module, env=rsaga.env(), show=TRUE)
{
    if (is.function(lib))
        lib = deparse(substitute(lib))
    
    if (substr(lib,1,6)=="rsaga.") {
        if (lib=="rsaga.fill.sinks") {
            warning("'rsaga.fill.sinks' uses three modules from the 'ta_preprocessor' library:\n",
                "   for 'method=\"planchon.darboux.2001\"': module 2\n",
                "   for 'method=\"wang.liu.2006\"': module 3\n",
                "   for 'method=\"xxl.wang.liu.2006\"': module 4\n",
                "using 'module=2'\n")
        }
        lib = switch(lib,
                rsaga.close.gaps          = list(lib="grid_tools", module=7),
                rsaga.esri.to.sgrd        = list(lib="io_grid", module=1),
                rsaga.sgrd.to.esri        = list(lib="io_grid", module=0),
                rsaga.parallel.processing = list(lib="ta_hydrology", module=0),
                rsaga.local.morphometry   = list(lib="ta_morphometry", module=0),
                rsaga.slope               = list(lib="ta_morphometry", module=0),
                rsaga.aspect              = list(lib="ta_morphometry", module=0),
                rsaga.curvature           = list(lib="ta_morphometry", module=0),
                rsaga.plan.curvature      = list(lib="ta_morphometry", module=0),
                rsaga.profile.curvature   = list(lib="ta_morphometry", module=0),
                rsaga.sink.route          = list(lib="ta_preprocessor", module=0),
                rsaga.sink.removal        = list(lib="ta_preprocessor", module=1),
                rsaga.fill.sinks          = list(lib="ta_preprocessor", module=2),
                rsaga.contour             = list(lib="shapes_grid", module=5),
                rsaga.hillshade           = list(lib="ta_lighting", module=0),
                #rsaga.solar.radiation     = list(lib="ta_lighting", module=2),
                #rsaga.insolation          = list(lib="ta_lighting", module=3),
                rsaga.filter.simple       = list(lib="grid_filter", module=0),
                rsaga.filter.gauss        = list(lib="grid_filter", module=1) )
        module = lib$module
        lib = lib$lib
    }

    res = NULL

    usage = rsaga.geoprocessor(lib, module, param = list(h=""), env = env,
        intern = TRUE, show.output.on.console = FALSE, flags = NULL,
        check.module.exists = FALSE, warn = -1)

    skip = 0
    while ((length(usage)>(1+skip)) & (substr(usage[1+skip],1,6)!="Usage:")) {
        if (substr(usage[1+skip],1,8) %in% 
                c("SAGA CMD","Copyrigh","library ","module n","________")) {
            skip = skip + 1
        } else {
            if (skip == 0) {
                usage = usage[ 2 : length(usage) ]
            } else {
                usage = usage[ c(1:skip, (skip+2):length(usage)) ]
            }
        }
    }
    if (length(usage) > 1) {
        res = usage[ 1 : (length(usage)-1) ]
        if (substr(res[length(res)],1,6)=="______") {
            res = c(res, "Usage description not available (interactive module?)")
            warning("usage description not available for module ",
                module, " in library ", lib, " (interactive module?)")
        }
    } else
        warning("usage description not available for module ",
            module, "\nin library ", lib, " (interactive module?)")
    if (show) {
        if (!is.null(res))
            cat(paste(res,collapse="\n"),"\n\n")
    }
    invisible(res)
}


#' HTML help on a SAGA module or library
#' 
#' This function opens SAGA's HTML documentation for the specified library or module. Works with SAGA GIS 2.1.0(+), for earlier versions a web page with the SAGA GIS wiki is displayed.
#'
#' @name rsaga.html.help
#' @param lib name of the SAGA library, or one of the \code{rsaga.} module functions such as \code{\link{rsaga.hillshade}}
#' @param module name or numeric identifier of SAGA module in library \code{lib}; \code{module=NULL} takes you to the main help page of the SAGA library \code{lib}
#' @param use.program.folder logical; if \code{TRUE} (the default), attempt to write SAGA GIS documentation to a \code{"help"} subfolder of \code{env$path}; the \code{"help"} folder is created if it doesn't exist. If \code{FALSE}, create SAGA GIS documentation files in this R session's temporary folder as obtained using \code{tempdir()}
#' @param env a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}
#' @param ... additional arguments to \code{\link{browseURL}}
#' @details Requires SAGA GIS 2.1.0(+), with earlier versions use \code{\link{rsaga.get.usage}}.
#' @examples
#' \dontrun{
#' # Requires SAGA GIS 2.1.0+:
#' rsaga.html.help("io_grid")
#' rsaga.html.help("io_grid",0)
#' rsaga.html.help("io_grid","Import ESRI Arc/Info Grid")
#' }
#' @seealso \code{\link{rsaga.get.usage}}, \code{\link{rsaga.geoprocessor}}, \code{\link{rsaga.env}}
#' @keywords utilities interface
#' @export
rsaga.html.help = function(lib, module=NULL, use.program.folder = TRUE, env=rsaga.env(), ...)
{
    if (env$version == "2.1.0" | env$version == "2.1.1" | env$version == "2.1.2" |
        env$version == "2.1.3" | env$version == "2.1.4" | env$version == "2.2.0" |
        env$version == "2.2.1" | env$version == "2.2.2" | env$version == "2.2.3") {
        # Convert character string module names to integer code, if possible:
        if (!is.null(module)) {
            if (is.character(module)) {
                modules = rsaga.get.lib.modules(lib, env=env)
                if (any(modules$name == module))
                    module = modules$code[ which(modules$name == module)[1] ]
            }
        }
    
        env$workspace = file.path(env$path, "help")

        # check if help file already exists:
        fnm = file.path(env$workspace,lib,lib)
        if (!is.null(module)) {
            if (is.numeric(module)) {
                fnm = paste(fnm, "_", formatC(module,width=2,flag="0"), sep="")
            }
        }
        fnm = paste(fnm, ".html", sep="")

        # Create help files in the SAGA GIS program folder if they don't yet exist,
        # or in this R session's temporary folder if need be:
        if(!file.exists(fnm)) {
            if (!file.exists(env$workspace)) {
                if (use.program.folder)
                    use.program.folder = dir.create(env$workspace)
                if (!use.program.folder) {
                    temp.workspace = file.path(tempdir(),"help")
                    if (!file.exists(temp.workspace)) dir.create(temp.workspace)
                    env$workspace = temp.workspace
                }
            }
        }
        # Updated path and file name (env$workspace may have changed):
        fnm = file.path(env$workspace,lib,lib)
        if (!is.null(module)) {
            if (is.numeric(module)) {
                fnm = paste(fnm, "_", formatC(module,width=2,flag="0"), sep="")
            }
        }
        fnm = paste(fnm, ".html", sep="")
        # Create help files in the SAGA GIS program folder if they don't yet exist:
        if (!file.exists(fnm)) {
            cat("Calling SAGA GIS to create documentation files in the program folder...\n")
            rsaga.geoprocessor(lib = NULL, module = NULL, prefix = "-d", env = env)
        }

        # Issue a warning if documentation file still can't be accessed:
        if (!file.exists(fnm)) {
            warning("Can't create or can't find suitable documentation file\n",
                    "for module ", module, " in library ", lib, ".\n",
                    "Possible reasons:\n",
                    "- No writing permission for folder ", env$workspace, "\n",
                    "- Module or library does not exist\n")
            return()
        }
        utils::browseURL( paste("file://",fnm,sep=""), ...)
    } else {
        warning("rsaga.html.help only available for SAGA GIS 2.1.0+.\n",
                "Redirecting you to the SAGA GIS wiki...")
        url = "http://sourceforge.net/apps/trac/saga-gis/wiki"
        utils::browseURL(url, ...)
    }
    return()
}


#' Generic R interface for SAGA modules
#'
#' This function is the workhorse of the R--SAGA interface: It calls the SAGA command line tool to run SAGA modules and pass arguments.
#'
#' @name rsaga.geoprocessor
#' @param lib Name of the SAGA library to be called (see Details).
#' @param module Number (\code{>=0}) or name of the module to called within the library \code{lib} (see Details).
#' @param param A list of named arguments to be passed to the SAGA module (see Examples).
#' @param show.output.on.console a logical (default: \code{TRUE}), indicates whether to capture the output of the command and show it on the R console (see \code{\link{system}}).
#' @param invisible a logical, indicates whether the command window  should be visible on the screen.
#' @param intern a logical, indicates whether to make the output of the command an R object
#' @param prefix optional character string: prefix such as \code{"-h"} used in the \code{saga_cmd} call; mostly for internal purposes; call \code{saga_cmd -h} from the command line for details; see also \code{flags}
#' @param flags optional character string indicating any command line flags; supported only by SAGA GIS 2.1.0 (and higher), quietly ignored otherwise: \code{"q"}: no progress report (the default for \code{show.output.on.console=TRUE}); \code{"r"}: no messages report; \code{"s"}: silent mode, i.e. no progress and no messages report  (the default for \code{show.output.on.console=FALSE}); other flag options probably not relevant within RSAGA
#' @param cores optional numeric argument, or \code{NA}: number of cores used by SAGA GIS; supported only by SAGA GIS 2.1.0 (and higher), ignored otherwise (with a warning); overwrites the \code{cores} setting specified in the \code{env} argument (see \code{\link{rsaga.env}}). Multicore-enabled SAGA GIS modules such as the one used by \code{\link{rsaga.pisr}} seem to run in multicore mode by default when this argument is not specified, therefore \code{cores} should only be specified to use a smaller number of cores than available on a machine.
#' @param env A SAGA geoprocessing environment, i.e. a list with information on the SAGA and SAGA modules paths and the name of the working directory in which to look for input and output files. (Defaults: see \code{\link{rsaga.env}}.)
#' @param display.command Display the DOS command line for executing the SAGA module (including all the arguments to be passed). Default: \code{FALSE}.
#' @param reduce.intern If \code{intern=TRUE}, reduce the text output of SAGA returned to R by eliminating redundant lines showing the progress of module execution etc. (default: \code{TRUE}).
#' @param check.module.exists logical (default: \code{TRUE}): call \code{\link{rsaga.module.exists}} to determine if the specified module can be called in the current SAGA installation
#' @param warn logical (default: \code{TRUE}): for internal purposes - can be used to suppress warning messages generated by failed SAGA_CMD calls; currently used by \code{\link{rsaga.get.lib.modules}} and related functions; see \code{\link{options}} argument \code{warn} for details
#' @param argsep character (default: \code{" "}; currently for internal use): defines the character symbol used as a separator between each argument name and argument value passed to \code{saga_cmd}. SAGA GIS 2.1.0 (RC1) seems to move toward \code{"="} as a separator, but \code{" "} still works and some modules (e.g. the used by \code{rsaga.pisr}) don't seem to work with \code{argsep="="}. Future releases of RSAGA may change the default \code{argsep} value and/or delete or ignore this argument and/or move it to \code{\link{rsaga.env}}.
#' @param ... Additional arguments to be passed to \code{\link[base]{system}}.
#' 
#' @details This workhorse function establishes the interface between the SAGA command line program and R by submitting a system call. This is a low-level function that may be used for directly accessing SAGA; specific functions such as \code{rsaga.hillshade} are intended to be more user-friendly interfaces to the most frequently used SAGA modules. These higher-level interfaces support default values for the arguments and perform some error checking; they should therefore be preferred if available.
#' 
#' A warning is issued if the RSAGA version is not one of 2.0.4-2.0.8 or 2.1.0-2.1.4
#'
#' @return The type of object returned depends on the \code{intern} argument passed to \code{\link{system}}.
#' 
#' If \code{intern=FALSE}, a numerical error/success code is returned, where a value of \code{0} corresponds to success and a non-zero value indicates an error. Note however that the function always returns a success value of \code{0} if \code{wait=FALSE}, i.e. if it does not wait for SAGA to finish.
#' 
#' If \code{intern=TRUE} (default), the console output of SAGA is returned as a character vector. This character vector lists the input file names and modules arguments, and gives a more or less detailed report of the function's progress. Redundant information can be cancelled out by setting \code{reduce.intern=TRUE}.
#'
#' @references Brenning, A., 2008. Statistical geocomputing combining R and
#'  SAGA: The example of landslide susceptibility analysis with
#'  generalized additive models. In J. Boehner, T. Blaschke and
#'  L. Montanarella (eds.), SAGA - Seconds Out (= Hamburger
#'  Beitraege zur Physischen Geographie und
#'  Landschaftsoekologie, vol. 19), p. 23-32.
#'
#' @author Alexander Brenning (R interface); Olaf Conrad and the SAGA development team (SAGA development)
#' @note Existing output files will be overwritten by SAGA without prompting!
#' 
#' If a terrain analysis function is not directly interfaced by one of the RSAGA functions, you might still find it in the growing set of SAGA libraries and modules. The names of all libraries available in your SAGA installation can be obtained using \code{\link{rsaga.get.libraries}} (or by checking the directory listing of the \code{modules} folder in the SAGA directory). The names and numeric codes of all available modules (globally or within a specific library) are retreived by \code{\link{rsaga.get.modules}}. Full-text search in library and module names is performed by \code{\link{rsaga.search.modules}}. For information on the usage of SAGA command line modules, see \code{\link{rsaga.get.usage}}, or the RSAGA interface function if available.
#' 
#' \code{display.command=TRUE} is mainly intended for debugging purposes to check if all arguments are passed correctly to SAGA CMD.
#' @seealso \code{\link{rsaga.env}}, \code{\link{rsaga.get.libraries}}, \code{\link{rsaga.get.modules}}, \code{\link{rsaga.search.modules}}, \code{\link{rsaga.get.usage}}; \code{\link{rsaga.esri.wrapper}} for a wrapper for ESRI ASCII/binary grids; \code{\link{rsaga.hillshade}} and other higher-level functions.
#' @examples
#' \dontrun{
#' rsaga.hillshade("dem","hillshade",exaggeration=2)
#' # using the RSAGA geoprocessor:
#' rsaga.geoprocessor("ta_lighting",0,list(ELEVATION="dem.sgrd",SHADE="hillshade",EXAGGERATION=2))
#' # equivalent DOS command line call:
#' # saga_cmd.exe ta_lighting 0 -ELEVATION dem.sgrd -SHADE hillshade -EXAGGERATION 2 
#' }
#' @keywords spatial interface
#' @export
rsaga.geoprocessor = function(
    lib, module = NULL, param = list(), 
    show.output.on.console = TRUE, invisible = TRUE, intern = TRUE,
    prefix = NULL, flags = ifelse(show.output.on.console,"q","s"), cores,
    env = rsaga.env(), display.command = FALSE, reduce.intern = TRUE,
    check.module.exists = TRUE, warn = options("warn")$warn,
    argsep = " ", ... )
{
# Issue warning if using SAGA GIS version that has not been tested with RSAGA:
    if (!is.null(env$version)) {
        if (!is.na(env$version)) {
            if (!any(c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8",
                       "2.1.0","2.1.1","2.1.2","2.1.3","2.1.4",
                       "2.2.0","2.2.1","2.2.2","2.2.3") == env$version))
                warning("This RSAGA version has been tested with SAGA GIS versions 2.0.4 - 2.2.2.\n",
                    "You seem to be using SAGA GIS ", env$version, ", which may cause problems due to\n",
                    "changes in names and definitions of SAGA module arguments, etc.", sep = "" )
        }
    }
    
    # Number of cores for multicore processing:
    if (!missing(cores)) env$cores = cores
    if (!is.na(env$cores) & (substr(env$version,1,4) == "2.0.")) {
        warning("'cores' argument not supported by SAGA GIS versions <2.1.0; ignoring 'cores' argument\n",
            "Use SAGA GIS 2.1.0+ for multicore geoprocessing.")
        env$cores = NA
    }

    if (is.na(env$version) | (substr(env$version,1,4) == "2.0.")) {
        if (argsep != " ")
            warning("To my knowledge, SAGA GIS 2.0.x only supports argsep=' '.\nUse different argsep values at own risk.")
    } else {
        if (!(argsep %in% c(" ","=")))
            warning("To my knowledge, SAGA GIS only supports argsep=' ' or (partially also) '='.\nUse different argsep values at own risk.")
    }
    
    # Change working directory:
    old.wd = getwd()
    on.exit(setwd(old.wd))
    setwd(env$workspace)

    # Set environment variables SAGA and SAGA_MLB:
    # (This might be redundant, but it probably won't hurt. Might also be version specific.)
    if ((Sys.info()["sysname"] != "Windows") | is.na(env$version) | (env$version == "2.0.4")) {
        old.saga = Sys.getenv("SAGA", unset = NA)
        old.saga.mlb = Sys.getenv("SAGA_MLB", unset = NA)
        on.exit(if (is.na(old.saga)) Sys.unsetenv("SAGA") else Sys.setenv(SAGA=old.saga), add = TRUE)
        on.exit(if (is.na(old.saga.mlb)) Sys.unsetenv("SAGA_MLB") else Sys.setenv(SAGA_MLB=old.saga.mlb), add=TRUE)
        Sys.setenv(SAGA=env$path, SAGA_MLB=env$modules)
    }

    # Core part of system call:    
    command = shQuote( paste( env$path, .Platform$file.sep, env$cmd, sep="" ) )

    # Prefix e.g. -h or --help for general help (this is currently used by rsaga.get.version)
    if (!is.null(prefix))
        command = paste( command, prefix, sep = " " )

    if (!(env$version %in% c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8","2.0.9"))) {
        if (!is.null(flags))
            command = paste( command, " -f=", flags, sep = "" )
    
        if (!is.null(lib) & !is.null(module) & (length(param)>0)) {
            if (!is.na(env$cores)) {
                command = paste( command, " --cores=", env$cores, sep="" )
            }
        }
    }
    
    # Library - in the case of unix systems (until 2.0.9), it must be preceded by 'lib' -
    # but not in the case of Mac OSX:
    ###add.lib = (Sys.info()["sysname"] != "Windows") & (Sys.info()["sysname"] != "Darwin")
    if (!is.null(lib)) {
        # From 2.1.0 on, UNIX-like systems do not have preceding 'lib' any more
        if (is.null(env$lib.prefix) |
            !(env$version %in% c("2.0.4", "2.0.5", "2.0.6",
                                 "2.0.7", "2.0.8", "2.0.9"))) env$lib.prefix = ""
        command = paste( command, " ", env$lib.prefix, lib, sep = "")
    }
    
    if (!is.null(lib) & !is.null(module)) {
        if (check.module.exists) {
            ex = rsaga.module.exists(lib, module, env=env)
            if (!ex) {
                cat("Module '", module, "' not found in SAGA library '", lib, "'.\n",
                     "Check if module name has changed (or is misspelled)?\n", sep = "")
                cat("The following (non-interactive) modules currently exist in this SAGA library:\n\n")
                print(rsaga.get.modules(lib, env=env, interactive = FALSE))
                cat("\n")
                stopifnot(rsaga.module.exists(lib,module, env=env))
            }
        }
    
        if (is.character(module)) module = shQuote(module)
        command = paste(command, module)
        if (length(param)>0) {
            # Logical arguments are treated in a special way with SAGA versions below 2.1.3:
            # They are simply omitted if their value is FALSE.
            if (!(env$version %in% c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8","2.0.9",
                                     "2.1.0","2.1.1","2.1.2"))) {
                i = 1
                while (i<=length(param)) {
                    if (is.logical(param[[i]])) {
                        if (!param[[i]]) {
                            param[[i]] = "false"
                            i = i - 1
                        } else param[[i]] = "true"
                    }
                    i = i + 1
                }
            } else {
                i = 1
                while (i<=length(param)) {
                    if (is.logical(param[[i]])) {
                        if (!param[[i]]) {
                            param[[i]] = NULL
                            i = i - 1
                        } else param[[i]] = ""
                    }
                    i = i + 1
                }
            }
            
            # Argument names:
            nm = names(param)
            # Argument values:
            val = as.character(unlist(param))
            
            # line added by Johan v.d.W.:
            # Put quotes around non-void argument values:
            val[ nchar(val) > 0 ] = shQuote( val[ nchar(val) > 0 ] )
            
            # Add saga_cmd arguments to the command line call:
            param = paste("-",nm, argsep, val,sep="",collapse=" ")
            command = paste(command, param)
        }
    }

    if (display.command) cat(command,"\n")

    if (Sys.info()["sysname"] == "Windows") {
        # Some rsaga core calls need to suppress warnings
        # related to non-zero exit codes of saga_cmd:
        oldwarn = options("warn")$warn
        on.exit(options(warn = oldwarn), add = TRUE)
        options(warn = warn)
        # Actual saga_cmd call:
        res = system( command, intern=intern,
            show.output.on.console=show.output.on.console, 
            invisible=invisible, ...)
        options(warn = oldwarn)
    } else {
        oldwarn = options("warn")$warn
        on.exit(options(warn = oldwarn), add = TRUE)
        options(warn = warn)
        res = system( command, intern=intern, ...)
        # 'show.output.on.console' and 'invisible' only work under Windows
        options(warn = oldwarn)
    }
    if (intern) {
        if (reduce.intern) {
            remove = grep("\r",res,fixed=TRUE)
            if (length(remove) > 0)
                res = res[ -remove ]
            remove = grep("^.*##.*##",res)
            if (length(remove) > 0)
                res = res[ -remove ]
            if (any(remove <- res=="go...")) res = res[!remove]
            if (any(remove <- res=="okay"))  res = res[!remove]
            if (any(remove <- substr(res,1,7)=="type -h")) res = res[!remove]
            if (any(remove <- substr(res,1,7)=="_______")) res = res[!remove]
        }
        if (show.output.on.console)
            cat(res,sep="\n")
    }
    if (intern) {
        invisible(res)
    } else   return(res)
}


#' Use RSAGA functions for ESRI grids
#'
#' This wrapper converts input grid files provided in ESRI binary (.flt) or ASCII (.asc) formats to SAGA's (version 2) grid format, calls the RSAGA geoprocessing function, and converts the output grids back to the ESRI grid format. Conversion can also be limited to either input or output grids.
#' @name rsaga.esri.wrapper
#' @param fun function: one of the RSAGA geoprocessing functions, such as \code{\link{rsaga.close.gaps}} or \code{\link{rsaga.hillshade}} etc.
#' @param in.esri logical: are input grids provided as ESRI grids (\code{in.esri=TRUE}) or as SAGA grids?
#' @param out.esri logical: should output grids be converted to ESRI grids?
#' @param env RSAGA environment as returned by \code{\link{rsaga.env}}
#' @param esri.workspace directory for the input and output ESRI ASCII/binary grids
#' @param format output file format, either \code{"ascii"} (default; equivalent: \code{format=1}) for ASCII grids or \code{"binary"} (equivalent: \code{0}) for binary ESRI grids (\code{.flt}).
#' @param georef character: \code{"corner"} (equivalent numeric code: \code{0}) or \code{"center"} (default; equivalent: \code{1}). Determines whether the georeference will be related to the center or corner of its extreme lower left grid cell.
#' @param prec number of digits when writing floating point values to ASCII grid files (only relevant if \code{out.esri=TRUE}).
#' @param esri.extension extension for input/output ESRI grids: defaults to \code{.asc} for \code{format="ascii"}, and to \code{.flt} for \code{format="binary"}
#' @param condensed.res logical: return only results of the RSAGA geoprocessing function \code{fun} (\code{condensed.res=TRUE}), or include the results of the import and export operations, i.e. the calls to \code{\link{rsaga.esri.to.sgrd}} and \code{\link{rsaga.sgrd.to.esri}}? (see Value)
#' @param clean.up logical: delete intermediate SAGA grid files?
#' @param intern \code{intern} argument to be passed to \code{\link{rsaga.geoprocessor}}; see Value
#' @param ... additional arguments for \code{fun}; NOTE: ESRI ASCII/float raster file names should NOT include the file extension (.asc, .flt); the file extension is defined by the \code{esri.extension} and \code{format} arguments!
#' @details ESRI ASCII/float raster file names should NOT include the file extension (.asc, .flt); the file extension is defined by the \code{esri.extension} and \code{format} arguments!
#' @return The object returned depends on the \code{condensed.res} arguments and the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}.
#'
#' If \code{condensed.res=TRUE} and \code{intern=FALSE}, a single numerical error code (0: success) is returned. If \code{condensed.res=TRUE} and \code{intern=TRUE} (default), a character vector with the module's console  output is returned (invisibly).
#'
#' If \code{condensed.res=FALSE} the result is a list with components \code{in.res}, \code{geoproc.res} and \code{out.res}. Each of these components is either an error code (for \code{intern=FALSE}) or  (for \code{intern=TRUE}) a character vector with the console output of the input (\code{\link{rsaga.esri.to.sgrd}}), the geoprocessing (\code{fun}), and the output conversion (\code{\link{rsaga.sgrd.to.esri}}) step, respectively. For \code{in.esri=FALSE} or \code{out.esri=FALSE}, the corresponding component is \code{NULL}.
#' @note Note that the intermediate grids as well as the output grids may overwrite existing files with the same file names without prompting the user. See example below.
#' @seealso \code{\link{rsaga.esri.to.sgrd}}, \code{\link{rsaga.sgrd.to.esri}}, \code{\link{rsaga.geoprocessor}}, \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' rsaga.esri.wrapper(rsaga.hillshade,in.dem="dem",out.grid="hshd",condensed.res=FALSE,intern=FALSE)
#' # if successful, returns list(in.res=0,geoproc.res=0,out.res=0),
#' # and writes hshd.asc; intermediate files dem.sgrd, dem.hgrd, dem.sdat,
#' # hshd.sgrd, hshd.hgrd, and hshd.sdat are deleted.
#' # hshd.asc is overwritten if it already existed.
#' }
#' @keywords spatial interface
#' @export
rsaga.esri.wrapper = function(fun, in.esri=TRUE, out.esri=TRUE, 
    env=rsaga.env(), esri.workspace=env$workspace,
    format="ascii", georef="corner", prec=5, esri.extension,
    condensed.res=TRUE, clean.up=TRUE, intern=TRUE, ...)
{
    in.res = NULL
    geoproc.res = NULL
    out.res = NULL
    format = match.arg.ext(format,choices=c("binary","ascii"),base=0,ignore.case=TRUE,numeric=TRUE)
    if (missing(esri.extension))
        esri.extension = c(".flt",".asc")[format+1]
    args = list(...)
    argnms = names(args)
    
    in.ok = TRUE
    if (in.esri) {
        wh = grep("^in\\.",names(args))
        if (length(wh)==0) {
            warning("'in.esri' is TRUE, but the geoprocessing function does not have an 'in.*' grid argument")
        } else {
            in.args = args[wh]
            in.res = rsaga.esri.to.sgrd(in.grids=set.file.extension(unlist(in.args),esri.extension),
                intern=intern, show.output.on.console=FALSE,
                out.sgrds=unlist(in.args), in.path=esri.workspace, env=env) # more args to geoproc
            if (!intern) in.ok = all(in.res==0)
        }
    }
    
    geoproc.ok = TRUE
    if (in.ok) {
        geoproc.res = fun(env=env,intern=intern,...)
        if (!intern) geoproc.ok = all(geoproc.res==0)
    }
    if (clean.up) {
        del.files = set.file.extension(in.args,"")
        del.files = unlist(lapply(as.list(del.files), function(x) paste(x,c("sgrd","hgrd","sdat","mgrd"),sep="")))
        unlink(del.files)
    }
    
    out.ok = TRUE
    if (out.esri & in.ok & geoproc.ok) {
        wh = grep("^out\\.",names(args))
        if (length(wh)==0) {
            warning("'out.esri' is TRUE, but the geoprocessing function does not have an 'out.*' grid argument")
        } else {
            out.args = args[wh]
            out.res = rsaga.sgrd.to.esri(in.sgrds=unlist(out.args),
                out.grids=set.file.extension(out.args,unlist(esri.extension)),
                out.path=esri.workspace, env=env, intern=intern, show.output.on.console=FALSE,
                format=format, georef=georef, prec=prec) # more args to geoproc
            if (!intern) out.ok = all(out.res==0)
            if (clean.up) {
                del.files = set.file.extension(out.args,"")
                del.files = unlist(lapply(as.list(del.files), function(x) paste(x,c("sgrd","hgrd","sdat","mgrd"),sep="")))
                unlink(del.files)
            }
        }
    }

    res = list( in.res=in.res, geoproc.res=geoproc.res, out.res=out.res )
    if (condensed.res) {
        if (intern) {
            res = geoproc.res
        } else   res = max(abs(unlist(res)))
    }
    if (intern) {
        invisible(res)
    } else  return( res )
}


#' Create a copy of a SAGA grid file
#' 
#' Creates a copy of a SAGA grid file, optionally overwriting the target file if it already exists. Intended mainly for internal use by RSAGA functions, currently in particular \code{\link{rsaga.inverse.distance}}.
#'
#' @param in.grid name of a SAGA GIS grid file; file extension can be omitted
#' @param out.grid name of a SAGA GIS grid file; file extension can be omitted
#' @param overwrite logical; if \code{TRUE} (the default), overwrite \code{out.grid} if it already exists; if \code{FALSE} and the \code{out.grid} already exists, copying will be skipped without causing an error.
#' @param env a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}
#' @note SAGA grid files consist of three (or more) individual files with file extensions \code{.mgrd}, \code{.sgrd} and \code{.sdat}. The files with these three file extensions are copied, any additional files (e.g. a history file) are ignored.
#' @keywords spatial interface
#' @export
rsaga.copy.sgrd = function(in.grid, out.grid, overwrite = TRUE, env = rsaga.env())
{
    in.grid = set.file.extension(in.grid,".sgrd")
    out.grid = set.file.extension(out.grid,".sgrd")
    stopifnot(in.grid != out.grid)
    old.wd = getwd()
    setwd(env$workspace)
    stopifnot(file.exists(in.grid))
    in.files = c(
        set.file.extension(in.grid,".mgrd"),
        set.file.extension(in.grid,".sgrd"),
        set.file.extension(in.grid,".sdat") )
    out.files = c(
        set.file.extension(out.grid,".mgrd"),
        set.file.extension(out.grid,".sgrd"),
        set.file.extension(out.grid,".sdat") )
    res = c()
    for (i in 1:length(in.files)) {
        res[i] = file.copy(in.files[i], out.files[i], overwrite = overwrite)
    }
    setwd(old.wd)
    return(all(res))
}
