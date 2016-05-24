# -----------------------------------------------------------
# autodetectPython
# ================
#' @title Autodetects the settings for Windows
#' 
#' @description Autodetects the settings needed to connect to
#'              the python dll file (\strong{only Windows}).
#' @param pythonExePath a string containing the path to "python.exe"
#'        (e.g. "C:\\Python34\\python.exe").
#' @return Returns a list containing the information necessary to
#'         connect to Python if a compatible Python, version was found,
#'         raises an error otherwise.
#' @examples
#' \dontrun{   
#'   autodetectPython()
#'   autodetectPython("C:\\Python27\\python.exe")
#' }
# -----------------------------------------------------------
autodetectPython <- function(pythonExePath=NULL){
    if (pyIsConnected()) stop("when Python is connected to R the function autodetectPython doesn't work!")
    pyHome <- NULL
    rArch <- if (grepl("i386", R.version$arch)) '32bit' else '64bit'
    # NOTE: guessPythonExePathWhere and guessPythonPathEnvironmentVariables
    #       do essentially the same same but the function where isn't 
    #       available pre Windows 2003 since many people still use 
    #       WinXp I just leave it as backup! 
    if (is.null(pythonExePath)){
        pyCandidates <- guessPythonExePathWhere()
        pyCandidates <- filterCandidatesByArch(pyCandidates, rArch)
        if (!is.null(pyCandidates)){
            pythonExePath <- pyCandidates[1]
        }else{
            pyCandidates <- guessPythonExePathEnvironmentVariables()
            pyCandidates <- filterCandidatesByArch(pyCandidates, rArch)
            if (!is.null(pyCandidates)) pythonExePath <- pyCandidates[1]
        }
    }
    msg <- paste(c('Python could not be found please provide the path to "python.exe"', 
                 '(e.g. "C:\\Python27\\python.exe") or set the appropriate path variable!')
                , collapse="\n\t")
    if (is.null(pythonExePath)) stop(msg)
 
    pyArch <- sprintf("%ibit", guessDllVersion(pythonExePath))
    if (pyArch != rArch) stop(sprintf("Python %s can't be connected with R %s!", pyArch, rArch))
    
    # an alternative would be to force the user to specify it
    # get major version
    pyMajorVersion <- as.integer(sysCallPython(pythonExePath, "str(sys.version_info.major)"))
    # get minor version
    pyMinorVersion <- as.integer(sysCallPython(pythonExePath, "str(sys.version_info.minor)"))                                      
                                               
    # For simplicity I will just assume PYTHONHOME is where python.exe
    # is located (Another approach would be to get python path and look
    # which folder is sub folder to most of the other folders. But if
    # if a user thinks it is a good idea to change the python folder
    # structure he should just specify the settings via pyConnectWinDll.
    # Also reading sys.prefix should normally work.)
    pyHome <- dirname(pythonExePath)
    
    dllName <- sprintf("python%i%i.dll", pyMajorVersion, pyMinorVersion)
    
    # it would be the easiest to get the dll from the win32api package
    # win32api.GetModuleFileName(sys.dllhandle) but one can not assume
    # that everyone has win32api installed
    # sysCallPython(pythonExePath, "win32api.GetModuleFileName(sys.dllhandle)", "import win32api;")
    #
    # For now I just look in the PYTHONHOME folder which yields to the 
    # following cases:
    # 1) portable python   
    # 2) if python is registered in the path windows will find it
    # 3) the user should register python in the path or provide the folder
    if (any(grepl(dllName, dir(pyHome), fixed=TRUE))){
        dllDir = pyHome
    }else{
        dllDir = NULL
    }
    list(pythonExePath=pythonExePath, 
         dllName=dllName, 
         dllDir=dllDir,
         majorVersion=pyMajorVersion,
         pythonHome=pyHome,
         arch=pyArch)
}

