# ------------------------------------------------------------------------------ 
#
#   PyRunString
#
# ------------------------------------------------------------------------------

#  -----------------------------------------------------------------------------
#  pyExecp
#  =======
#' @title Executes a single line of Python code from within R
#' @description The function pyExecp is designed to execute a single line of 
#'              Python code from within R. Thereby pyExecp tries to emulate
#'              the natural interactive Python terminal behavior.
#' @param code a string of Python code to be executed in Python.
#' @details The name pyExecp is short for python execute and print.
#'          As the name is indicating the most visual difference between pyExec 
#'          and pyExecp lies in the printing behavior. For example, executing 
#'          \code{pyExecp('"Hello " + "R!"')} would show \code{'Hello R!'} in the
#'          R terminal, executing \code{pyExec('"Hello " + "R!"')} wouldn't show 
#'          anything. Internally pyExec uses PyRun_SimpleString and pyExecp uses 
#'          PyRun_String with the flag Py_single_input, therefore pyExecp 
#'          can be used to simulate an interactive Python interpreter behavior.
#'
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' pyExecp('"Hello" + " " + "R!"')
#  -----------------------------------------------------------------------------   
pyExecp <- function(code){
    ret <- -1
    check_string(code)
    if (nchar(code) > 0){
        if ( pyConnectionCheck() ) return(invisible(NULL))
        if (pyOptions("winPython364")){
            ret <- try(.Call("py_run_string_single_input", code), silent = TRUE)
            cat(pyGetSimple("__getStdout()")) ## print stdout
            if (ret == -1){
                msg <- makeErrorMsg()
                if (!is.null(msg)) stop(msg)
            }
        }else{
            ret <- .Call("py_run_string_single_input", code)
            if (ret == -1) stop("An error has occured while executing Python code.",
                                " See traceback above.")
        }
    }
    return( invisible(ret) )
}

#  -----------------------------------------------------------------------------
#  pyExec
#  ======
#' @title Executes multiple lines of Python code from within R
#'
#' @description The function pyExec allows to execute multiple lines of python 
#'              code from within R. 
#' @param code a string of Python code to be executed in Python.
#' @details Since pyExec can execute multiple lines, it is the obvious choice for defining
#'          Python functions or running small scripts where no return value is needed.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' pyExec('
#' print("The following line will not appear in the R terminal!")
#' "Hello" + " " + "R!"
#' print("NOTE: pyExecp would also show the line above!")
#' print("The following line will appear in the R terminal!")
#' print("Hello" + " " + "R!")
#' ')
#  -----------------------------------------------------------------------------  
pyExec <- function(code){
    ret <- -1
    check_string(code)
    if (nchar(code) > 0){
        if ( pyConnectionCheck() ) return(invisible(NULL))
        if (pyOptions("winPython364")){
            ret <- try(.Call("py_run_simple_string", code), silent = TRUE)
            cat(pyGetSimple("__getStdout()")) ## print stdout
            if (ret == -1){
                msg <- makeErrorMsg()
                if (!is.null(msg)) stop(msg)
            }
        }else{
            ret <- .Call("py_run_simple_string", code)
            if (ret == -1) stop("An error has occured while executing Python code.",
                                " See traceback above.")
        }
    }
    return( invisible( ret ) )
}

## TODO: (not important)
## Was ich eigentlich haben moechte ist ein zusaetzlicher
##       parameter conversionErrors der wie bei Python Unicode entweder
##       die Werte "strict", "replace" und "ignore" annehmen kann.
##       "strict" raises an error
##       "replace" replaces the object by the string representation 
##                 and raises a warning
##       "ignore" replaces the object by the string representation no
##                warning or error gets
#  -----------------------------------------------------------------------------
#  pyExecg
#  =======
#' @title Executes multiple lines of python code and gets the output
#'
#' @description The function pyExecg is designed to execute multiple lines of 
#'              Python code and returns the thereby generated variables to R.
#' @param code a string of Python code to be executed in Python.
#' @param returnValues a character vector containing the names of the variables,
#'        which should be returned to R.
#' @param autoTypecast a an optional logical value, default is TRUE, specifying
#'        if the return values should be automatically typecasted if possible.
#' @param returnToR an optional logical, default is TRUE, specifying if the 
#'        generated variables should be returned to R.
#' @param mergeNamespaces an optional logical, default is FALSE, specifying if 
#'        the internally generated temporary namespace should be merged with the 
#'        name space __main__. See \bold{Details}.
#' @param override an optional logical value, default is FALSE, specifying how to 
#'        merge the temporary namespace with the __main__ namespace.
#' @param simplify an optional logical, if TRUE (default) R converts Python lists 
#'        into R vectors whenever possible, else it translates Python lists 
#'        always to R lists.
#' @return Returns a list containing all the variables of the __main__ namespace.
#' @details The function pyExecg executes the code in a temporary namespace, 
#'          after the execution every variable from the namespace is returned 
#'          to R. If the mergeNamespaces is set to TRUE the temporary namespace
#'          gets merged with the (global) namespace __main__.
#'          The logical variable override is used to control, if already
#'          existing variables in the namespace __main__ should be overridden,
#'          when a variable with the same name get's assigned to the
#'          temporary namespace. If a python object can't be converted to an
#'          R object it is assigned to the Python dictionary __R__.namespace
#'          and the type, id and an indicator if the object is a callable are 
#'          returned.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' if ( pyIsConnected() ){
#' # 1. assigns x to the global namespace
#' pyExec("x=4")
#' # 2. assigns y to the temp namespace
#' pyExecg("y=4", simplify=TRUE)
#' # 3. assign again to the temp namespace
#' pyExecg("
#' y=[i for i in range(1,4)]
#' x=[i for i in range(3,9)]
#' z=[i**2 for i in range(1,9)]
#' ", returnValues=c("x", "z"), simplify=TRUE)
#' # 4. assign x to the temp namespace, x gets returned as vector
#' pyExecg("x=[i for i in range(0,5)]", simplify=TRUE)
#' # 5. assign x to the temp namespace, x gets returned as list
#' pyExecg("x=[i for i in range(0,5)]", simplify=FALSE)
#' # 6. x is still 4 since except assignment 1 all other assignments
#' #    took place in the temp namespace
#' pyPrint("x")
#' # 7. note y has never been assigned to the main namespace
#' "y" %in% pyDir()
#' # 8. since mergeNamespaces is TRUE PythonInR will try
#' #     to assign x to the main namespace but since override is
#' #     by default FALSE and x already exists in the main namespace
#' #     x will not be changed
#' pyExecg("x=10", simplify=TRUE, mergeNamespaces=TRUE)
#' # 9. there is no y in the main namespace therefore it can be assigned
#' pyExecg("y=10", simplify=TRUE, mergeNamespaces=TRUE)
#' pyPrint("x") # NOTE: x is still unchanged!
#' pyPrint("y") # NOTE: a value has been assigned to y!
#' # 10. since override is now TRUE the value of x will be changed in the
#' #      main namespace
#' pyExecg("x=10", simplify=TRUE, mergeNamespaces=TRUE, override=TRUE)
#' pyPrint("x") # NOTE: x is changed now!
#' # 11. get an object which can't be typecast to an R object
#' #     pyExecg does not transform these objects automatically
#' pyExec("import os")
#' z <- pyExecg("x = os")
#' os <- PythonInR:::pyTransformReturn(z[[1]])
#' os$getcwd()
#' }
#  -----------------------------------------------------------------------------
pyExecg <- function(code, returnValues=character(), autoTypecast=TRUE, returnToR=TRUE, 
                    mergeNamespaces=FALSE, override=FALSE, simplify=TRUE){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    if (!is.character(returnValues)){
        stop('argument returnValue must be a character vector')
    }
    if (length(returnValues) > 0){
        rv <- sprintf("(%s)", paste(shQuote(returnValues), collapse=", "))
        del <- sprintf("
for __d00d__ in dir():
  if (__d00d__ not in %s):
    try:
      del(locals()[__d00d__])
    except: 
      pass

try:
  del(__d00d__)
except:
  pass
", rv)
        code <- sprintf("%s\n%s", code, del)
    }

    returnToR <- if (returnToR) 2L else 0L
    
    if (pyOptions("winPython364")){
        ret_val <- try(.Call("PythonInR_Run_String", code, 257L, autoTypecast,
                             mergeNamespaces, override, returnToR, 
                             simplify), silent = TRUE)
        cat(pyGetSimple("__getStdout()")) ## print stdout
        msg <- makeErrorMsg()
        if (!is.null(msg)) stop(msg)
    }else{
        ret_val <- .Call("PythonInR_Run_String", code, 257L, autoTypecast,
                         mergeNamespaces, override, returnToR, simplify)
    }
    # NOTE: the flag returnToR makes also a difference at the c level
    #       if it is FALSE only NULL get's returned!
    if (returnToR) return(ret_val)
    invisible(ret_val)
}

# An intern smaller version else I get a endless recursion when 
# I add printStoutErr also to pyExecg
pyExecgIntern <- function(code, autoTypecast=TRUE, mergeNamespaces=FALSE,
                          override=TRUE, simplify=TRUE){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    ret_val <- try(.Call("PythonInR_Run_String", code, 257L, autoTypecast,
                     mergeNamespaces, override, 2L, simplify), silent = TRUE)
    return(ret_val)
}

#  -----------------------------------------------------------------------------
#  pyExecfile
#  ==========
#' @title Executes Python source file from within R
#'
#' @description The function pyExecfile calls the Python function \code{execfile}. which 
#'              is the Python equivalent to the function \code{source} provided in R.
#' @param filename a character string giving the name or full path of the file to
#'                 be executed.
#' @details The function execfile is kind of the source of Python. Since it got omitted 
#'          in Python 3 a replacement gets assigned following the suggestions from \cr
#'          \url{http://www.diveintopython3.net/porting-code-to-python-3-with-2to3.html}.
#' @examples
#' \dontrun{
#' pyExecfile("myPythonScript.py")
#' }
#  ----------------------------------------------------------------------------- 
pyExecfile <- function(filename){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    check_string(filename)
    x <- pyExecp(sprintf("execfile('%s')", as.character(filename)[1]))
    return(invisible(x))
}

#  -----------------------------------------------------------------------------
#  pyPrint
#  =======
#' @title Convenience function to print a given Python object to the R terminal
#'
#' @description Prints a given Python variable to the R terminal.
#' @param objName a character string to be evaluated in Python and printed 
#'                to the R terminal.
#' @details Internally it uses a combination of pyExec and print. Please note that
#'          the result of pyExecp("x") and pyPrint("x") often will be different, 
#'          since pyExecp("x") is equivalent to typing x into the Python terminal
#'          whereas pyPrint("x") is equivalent to typing print(x) into the Python 
#'          terminal.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' pyPrint("'Hello ' + 'R!'")
#' pyPrint("sys.version")
#  ----------------------------------------------------------------------------- 
pyPrint <- function(objName){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    cmd <- sprintf("print(%s)", objName)
    return(invisible(pyExec(cmd)))
}
