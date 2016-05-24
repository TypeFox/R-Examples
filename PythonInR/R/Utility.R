# ------------------------------------------------------------------------------ 
#
#   Utility.r
#
# ------------------------------------------------------------------------------

pyIsCallable <- function(x){
    ## TODO: Test for 3 <= python < 3.2
    ## NOTE: I don't use try execpt since I wan't to see
    ##       the error if one occurs
    pyGet(sprintf('callable(%s)', x))
}

## if it exists it should have a type
pyVariableExists <- function(key){
cmd <- '
try:
    x = type(%s)
    x=True
except:
    x=False
'
    pyExecg(sprintf(cmd, key))[['x']]
}

## a small function to check if a variable name or a expression
isVariableName <- function(x){
    x = gsub("^[a-zA-Z\\_][a-zA-Z0-9\\_\\.]*", "", x, perl=TRUE)
    nchar(x) == 0
}

isBasic <- function(x){
    if (!is.null(dim(x))) return(FALSE)
    if (length(x) > 1) return(FALSE)
    if (is.list(x)) return(FALSE)
    TRUE
}

## variable
## variable name
## sollType
checkType <- function(penv, parentInfo, ...){
    x <- list(...)
    nam <- names(x)
    type <- as.character(x)
    b <- nam %in% ls(penv)
    if ( !all(b) ){
        varNames <- paste(nam[!b], collapse="', '")
        stop("the following variable names could not be found '",
             varNames, "'", call.=FALSE)
    }
    errMsg <- "in %s: argument '%s' must be of type 's'"
    for (i in 1:length(nam)){
        if (! inherits(penv[[nam[i]]], type[i]) ){
            stop(sprintf(errMsg, parentInfo, nam[i], type[i]), call.=FALSE)
        }
    }
    return(NULL)
}

# check_string 
#    checks the provided
## TODO: (change this to something nicer)
##       the function is not conistent any more
check_string <- function(x, minlen=1){
    vname <- deparse(substitute(x))
    errMes <- sprintf('argument "%s" must be a character vector of length 1', vname)
    if ( typeof(x) !=  "character" | length(x) != 1) stop(errMes)
    errMes <- sprintf('argument "%s" must have at least %i character', vname, minlen)
    if ( nchar(x) < minlen) stop(errMes)
}

# printStoutErr
#     Is a small work around wich fixes the issue with Python 3 64-bit and MinGW
makeErrorMsg <- function(){
    err = pyGetSimple('__getStderr()')
    if ( !is.null(err) ){
        if (nchar(err) > 0){
            # compile a error message and raise an error
            return(paste(c("", sprintf("   %s", unlist(strsplit(err, '\n')))), collapse="\n"))
        }
    }
    return(NULL)
}

# sysCallPython("C:\\Python27\\python.exe", "str(sys.version_info.major)", "import win32api;")
sysCallPython <- function(python, cmd, import=""){
    pyArgs <- "-c" 
    baseCmd <- "%simport sys;import os;sys.stdout.write(%s)"
    cmd <- sprintf(baseCmd, import, cmd)
    cmd <- paste(python, pyArgs, shQuote(cmd))
    returnValue <-
    tryCatch({
        system(cmd, intern=TRUE, ignore.stderr = TRUE)       
    }, warning = function(w) {
        NULL
    }, error = function(e) {
        NULL
    })
    returnValue
}

guessDllVersion <- function(dllPath){
    f <- file(dllPath, "rb")
    if (readChar(f, 2) != "MZ") return(-2)
    seek(f, 60, rw="rb")
    b = readBin(f, "raw", 4)
    header_offset <- unpack("V", b)[[1]]
    seek(f, header_offset+4, rw="rb")
    b = readBin(f, "raw", 2)
    bv <- unpack("V", b)[[1]]
    bit <- -1
    if (bv %in% c(332)){
        bit <- 32
    }else if (bv %in% c(512, 34404)){
        bit <- 64
    }
    close(f)
    return(bit)
}

filterCandidatesByArch <- function(pyCandidates, rArch){
    if ( length(pyCandidates) == 0 ) return(NULL)
    pyArchs <- paste(sapply(pyCandidates, guessDllVersion), "bit", sep="")
    pyCandidates <- pyCandidates[rArch == pyArchs]
    if (length(pyCandidates) == 0) return(NULL)
    pyCandidates
}

# Guess the path to python.exe utilizing the Windows batch 
# function where
# Returns: A character vector containing candidates, on success
#          NULL otherwise.
guessPythonExePathWhere <- function(){
    tryCatch({
        system("where python.exe", intern = TRUE, 
               ignore.stdout = FALSE, ignore.stderr = TRUE)
    }, warning = function(w) {
        NULL
    }, error = function(e) {
        NULL
    })
}

# Guess the path to python.exe utilizing the locations in 
# the environment variable path.
# Returns: A character vector containing candidates, on success
#          NULL otherwise.
guessPythonExePathEnvironmentVariables <- function(){
    pythonExePaths <- NULL
    path <- unlist(strsplit(Sys.getenv("PATH"), ";", fixed=TRUE))
    pyCandidates <- grep("python", path, value=TRUE, ignore.case=TRUE)
	if (length(pyCandidates) == 0) return(NULL)
    fun <- function(x) any(grepl("^python.exe", dir(x), ignore.case=TRUE))
    b <- sapply(pyCandidates, fun)
    pyCandidates <- pyCandidates[b]
    if (sum(b) > 0){
        # some one could get the idea to rename python.exe to Python.exe
        fun <- function(x) grep("^python.exe", dir(x), ignore.case=TRUE, value=TRUE)
        pyExeNames <- sapply(pyCandidates, fun)
        pythonExePaths <- normalizePath(file.path(pyCandidates, pyExeNames))
    }
    pythonExePaths
}

askPythonForPythonExePath <- function(){
    pythonPath <- 
    tryCatch({system(quote('python -c "import sys;import os;sys.stdout.write(os.path.dirname(sys.executable))"'), intern=TRUE)
    }, warning = function(w) {
        NULL
    }, error = function(e) {
        NULL
    })
    if (is.null(pythonPath)) return(NULL)
    fun <- function(x) grep("^python.exe", dir(x), ignore.case=TRUE, value=TRUE)
    pythonExeName <- fun(pythonPath)
    if (length(pythonExeName) > 0){
        return(normalizePath(file.path(pythonPath, pythonExeName)))
    }
    NULL
}

