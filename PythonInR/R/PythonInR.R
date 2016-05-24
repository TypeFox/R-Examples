# ------------------------------------------------------------------------------
#
#   PythonInR provides functions to interact with Python from within R.
#
# ------------------------------------------------------------------------------
registerLogCatcher <- function(){
code <- '
import io
sys_stdout = io.StringIO()
sys_stderr = io.StringIO()
sys_stdin = io.StringIO()
sys.stdout = sys_stdout
sys.stderr = sys_stderr
sys.stdin = sys_stdin

def __getStdout():
    x = sys.stdout.getvalue()
    sys_stdout.truncate(0)
    sys_stdout.seek(0)
    return(x)

def __getStderr():
    x = sys.stderr.getvalue()
    sys_stderr.truncate(0)
    sys_stderr.seek(0)
    return(x)

'
ret <- .Call( "py_run_simple_string",  code)
if ( ret == -1 ) warning("stdout redirection was not successful")
return( invisible( ret ) )
}

pyConnectWinDll <- function(dllName, dllDir, majorVersion, 
                            pythonHome, pyArch, useCstdout=NULL){
    if(pyIsConnected()){
        cat("R is already connected to Python!\n")
        return(NULL)
    }
    if (is.null(pyArch)) stop("couldn't detect Python architecture!")
    useAlteredSearchPath <- if (is.null(dllDir)) FALSE else TRUE
    if ((!is.null(dllName)) & (!is.null(dllDir))){
        dllPath <- file.path(dllDir, dllName)
        dllVersion <- guessDllVersion(dllPath)
        if (grepl("i386", R.version$arch)){
            if (dllVersion != 32) stop("64 bit dll (Python) can not be linked to 32 bit R!")
        }else{
            if (dllVersion != 64) stop("32 bit dll (Python) can not be linked to 64 bit R!")
        }
    }
    if (is.null(majorVersion)){
        majorVersion = as.integer(regmatches(dllName, regexpr("[0-9]", dllName)))
    }
    if ( is.null(useCstdout) ){
        useCstdout <- if ( (majorVersion==3) & (pyArch=='64bit') ) FALSE else TRUE
    }
    #cat("useCstdout: ", useCstdout, "\n")
    if (!is.null(pythonHome)){
        Sys.setenv(PYTHONHOME=pythonHome)
    }    
    .Call( "py_set_major_version", majorVersion)
    .Call( "py_connect", dllName, dllDir, as.integer(useAlteredSearchPath) )
    .Call( "py_get_process_addresses" )
    .Call( "py_set_program_name", "PythonInR" )
    if(useCstdout){
        .Call( "py_import_append_logCatcher" ) # has to be before py_initialize
    }
    .Call( "py_initialize", 1 )
    .Call( "py_init_py_values" )
    if(useCstdout){
        .Call( "py_init_redirect_stderrout" )
    }
    # import define a alternative to execfile as sugested at various sources
    # http://www.diveintopython3.net/porting-code-to-python-3-with-2to3.html
    if (majorVersion >= 3){
        if (!useCstdout){
            registerLogCatcher()
            pyOptions("winPython364", TRUE)
        }
        pyExec("def execfile(filename):\n    exec(compile(open(filename, 'rb').read(), filename, 'exec'), globals())");
    }
}

pyConnectStatic <- function(){
    .Call( "py_connect", 1 )
}

pyCranConnect <- function(){
    silent <- tryCatch({pyConnect()}, 
        error=function(e){
            cat(e$message)
            ## add some additional debug printing
            cat("print path:\n")
            print(Sys.getenv("PATH"))
        })
##    if ( Sys.info()['sysname'] == "Windows" ){
##        if (grepl("i386", R.version$arch)){
##            pyConnect()
##        }else{
##            # NOTE: Python 64bit is not available on the windows version of cran.
##            #       I just never connect to Python therefore the tests can run 
##            #       on the 32bit version as planed and I have not much to change 
##            #       if Python 64bit gets available.
##            #pyConnect()
##        }
##    }else{
##        # Since under Linux the default is static linkage which needs
##        # which needs no parameters.
##        pyConnect()
##    }
    invisible(NULL)
}

#  -----------------------------------------------------------------------------
#  pyConnect
#  =========
#' @title connects R to Python
#'
#' @description Connects R to Python. 
#'              \strong{(The parameters are only needed for the Windows version!)}
#' @param pythonExePath a character containing the path to "python.exe" 
#'                      (e.g. "C:\\Python27\\python.exe")
#' @param dllDir an optional character giving the path to the dll file. 
#'               Since the dll file is normally in a system folder or in the same 
#'               location as python.exe, this parameter is \bold{almost never needed}!
#' @param pythonHome an optional character giving the path to PYTHONHOME.
#'                   On Windows by default PYTHONHOME is the folder where python.exe 
#'                   is located, therefore this parameter is \bold{normally not needed}.
#' @details There is a different behavior for the static (Linux default)
#'          and the explicit linked (Windows default) version. Where as the
#'          static linked version automatically connects, when the package get's loaded,
#'          the explicitly linked version needs to be connected manually. 
#'          More information can be found at the README file or at 
#'          \url{http://pythoninr.bitbucket.org/}.
#' @note See the \href{https://cran.r-project.org/web/packages/PythonInR/README.html}{README} for more information about the Windows setup.
#' @examples
#' \dontrun{
#' ## Linux examples
#' pyConnect() # is done by default when the package is loaded
#'
#' ## Windows examples
#' pyConnect() ## will try to detect a suitable python version 
#'             ## from the PATH given in the environment variables
#' pyConnect("C:\\Python27\\python.exe")
#'
#' ## One can also explicitly set the parameters for the connection.
#' PythonInR:::pyConnectWinDll(dllName="python27.dll", dllDir=NULL,
#'                             majorVersion=2, pythonHome="C:\\Python27", 
#'                             pyArch="32bit")
#' }
#  -----------------------------------------------------------------------------
pyConnect <- function(pythonExePath=NULL, dllDir=NULL, pythonHome=NULL){
    if(pyIsConnected()){
        cat("R is already connected to Python!\n")
    }else{
        if (.Call( "isDllVersion")){
            py <- autodetectPython(pythonExePath)
            dllName <- py[['dllName']]
            if (is.null(dllDir)){
                dllDir <- py[['dllDir']]
            }else{
                if (!any(grepl(dllName, dir(dllDir)))){
                    stop(sprintf('"%s" could not be found at at the specified dllDir:\n\t"%s"!', dllName, dllDir))
                }
            }
            majorVersion <- py[['majorVersion']]
            pythonHome <- py[['pythonHome']]
            pyArch <- py[['arch']]
            silent <- pyConnectWinDll(dllName, dllDir, majorVersion, pythonHome, pyArch)
        }else{
            silent <- pyConnectStatic()
        }
        pyImportPythonInR()

        packageStartupMessage(sprintf("\nInitialize Python Version %s\n", pyVersion()))            

    }
    invisible(NULL)
}

#  -----------------------------------------------------------------------------
#  pyIsConnected
#  =============
#' @title checks if R is connected to Python
#'
#' @description Checks if R is connected to Python.
#' @return Returns TRUE if R is connected to Python, FALSE otherwise.
#' @examples
#' pyIsConnected()
#  -----------------------------------------------------------------------------
pyIsConnected <- function() as.logical(.Call( "py_is_connected" ))

# Prints some information about the Python version R is connected to.
# @examples
# \dontshow{PythonInR:::pyCranConnect()}
# pyInfo()
# Was once more use full but to be compatible with Python 3 it is
# now kind of useless.
pyInfo <- function(){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    .Call( "py_get_info" )
}

#  -----------------------------------------------------------------------------
#  pyVersion
#  =========
#' @title is a convenience function to get sys.version from Python
#'
#' @description A convenience function to get sys.version.
#' @return Returns a string containing the Python version and some compiler information.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' pyVersion()
#  -----------------------------------------------------------------------------
pyVersion <- function(){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    pyGet("sys.version")
}

#  -----------------------------------------------------------------------------
#  pyExit
#  ======
#' @title closes the connection to Python
#'
#' @description Closes the connection from R to Python.
#' @examples
#' \dontrun{
#' pyExit()
#' }
#  -----------------------------------------------------------------------------
pyExit <- function(){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    .Call( "py_close" )
}

# Checks if Python is connected to Python and prints a warning
# and returns TRUE if Python is not connected.
# This function only exists that I have to change the warning
# message only once if I would like to change it.
pyConnectionCheck <- function(){
    if ( !pyIsConnected() ){
        warning("R isn't connected to Python!\n")
        return(TRUE)
    }
    return(FALSE)
}
