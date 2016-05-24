## ----------------------------------------------------------------------------- 
##
##   PythonObjects
##
##  
## -----------------------------------------------------------------------------

pyObjectFinalize <- function(self){
    pyExec(sprintf("del(%s)", self$py.variableName))
}

callFun <- '
function(...){
  x <- list(...)
  i <- if ( !is.null(names(x)) ) (nchar(names(x)) > 0) else rep(FALSE, length(x))
  xargs <- if ( sum(!i) > 0 ) x[!i] else NULL
  xkwargs <- if ( sum(i) > 0 ) x[i] else NULL
  pyCall("%s", args=xargs, kwargs=xkwargs)
}
'

activeFun <- '
function(value){
    if (missing(value)){
        return(pyGet0("%s.%s"))
    }else{
        pySet("%s", value, "%s")
    } 
}
'

activeFun0 <- '
function(value){
    if (missing(value)){
        return(pyGet0("%s"))
    }else{
        pySet("%s", value)
    } 
}
'

## In Python try except is faster than if.
pyGetName <- function(x){
    pyExecg(sprintf('
try:
    x = %s.__name__
except:
    x = None
', x))[['x']]
}

#  ---------------------------------------------------------
#  pyObject
#  ========
#' @title Creates a virtual Python object
#'
#' @description The function pyObject creates a virtual Python object 
#'              of type PythonInR_Object.
#' @param key a character string giving the name of the Python object.
#' @param regFinalizer a logical indicating if a finalizer should be
#'                     be registered, the default value is TRUE.
#' @details Every PythonInR_Object has the following members:
#' \itemize{
#'   \item \strong{py.variableName} the variable name used in Python.
#'   \item \strong{py.objectName} the name of the Python object (obtained by x.__name__) 
#'          or NULL.
#'   \item \strong{py.type} the type of the Python object.
#'   \item \strong{py.del} a function to delete the Python object.
#'   \item \strong{print} for more information see R6 classes.
#'   \item \strong{initialize} for more information see R6 classes.
#' }
#' 
#' The other members of PythonInR_Object's are generated dynamically
#' based on the provided Python object. The R function \strong{ls} can be used
#' to view the members of a PythonInR_Object object.
#'
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' if ( pyIsConnected() ){
#' pyExec("import os")
#' os <- pyObject("os", regFinalizer = FALSE)
#' ls(os)
#' ## To show again the difference between pyGet and pyGet0.
#' os1 <- pyGet0("os") ## has no finalizer
#' os2 <- pyGet("os")  ## has a finalizer
#' os$py.variableName
#' os1$py.variableName
#' os2$py.variableName
#' }
#  ---------------------------------------------------------
pyObject <- function(key, regFinalizer = TRUE){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    check_string(key)

    objectName <- pyGetName(key)
    type <- pyType(key)

    pyMethods <- list()
    pyActive <- list()
    
    pydir <- pyDir(key)
    for (o in pydir){
        po <- paste(c(key, o), collapse=".")
        if (pyIsCallable(po)){
            cfun <- sprintf(callFun, po)
            pyMethods[[o]] <- eval(parse(text=cfun))
        }else{
            afun <- sprintf(activeFun, key, o, o, key)
            pyActive[[o]] <- eval(parse(text=afun))
        }
    }

    ## Choose names with a '.' since a point would violate the python
    ## name convention! This leaves me to take care of initialize and
    ## print where I can't chane the name. Therefore if a object 
    ## has a member with the name print it is renamed to py.print
    ## and initialize to py.initialize
    for (n in c("print", "initialize")){
        names(pyMethods)[names(pyMethods) == n] <- sprintf("py.%s", n)
        names(pyActive)[names(pyActive) == n] <- sprintf("py.%s", n)
    }

    if ( (!is.null(objectName)) & (!is.null(type)) ){
        className <- sprintf("%s.%s", type, objectName)
    }else if (is.null(objectName)){
        className <- type
    }else if (is.null(type)){ # should never happen since everything should have a type
        className <- objectName
    }else{
        className <- "?"
    }

    if (regFinalizer){
        pyobject <- R6Class(className,
                    portable = TRUE,
                    inherit = PythonInR_Object,
                    public = pyMethods,
                    active = pyActive)
    }else{
        pyobject <- R6Class(className,
                    portable = TRUE,
                    inherit = PythonInR_ObjectNoFinalizer,
                    public = pyMethods,
                    active = pyActive)
        class(pyobject) <- class(pyobject)[-2]
    }

    pyobject$new(key, objectName, type)
}

PythonInR_Object <- R6Class(
    "PythonInR_Object",
    public=list(
        portable=TRUE,
        py.variableName=NA,
        py.objectName="",
        py.type="",
        py.del = function(){
            pyExec(sprintf("del(%s)", self$py.variableName))
        },
        initialize = function(variableName, objectName, type) {
            if (!missing(variableName)) self$py.variableName <- variableName
            if (!missing(objectName)) self$py.objectName <- objectName
            if (!missing(type)) self$py.type <- type
            reg.finalizer(self, pyObjectFinalize, onexit = TRUE)
        },
        # #print = function(){pyExecp(self$py.variableName)}
        ## This should better handle unicode.
        print = function() pyPrint(self$py.variableName)
        ))

PythonInR_ObjectNoFinalizer <-
    R6Class("PythonInR_Object",
            portable = TRUE,
            inherit = PythonInR_Object,
            public = list(
                initialize = function(variableName, objectName, type) {
                    if (!missing(variableName)) self$py.variableName <- variableName
                    if (!missing(objectName)) self$py.objectName <- objectName
                    if (!missing(type)) self$py.type <- type
                }
            ))

#  ---------------------------------------------------------
#  pyFunction
#  ==========
#' @title creates a virtual Python function
#'
#' @description The function pyFunction creates a new object of type 
#'              pyFunction based on a given key.
#' @param key a string specifying the name of a Python method/function.
#' @details The function pyFunction makes it easy to create interfaces 
#'          to Python functions.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' if ( pyIsConnected() ){
#' pySum <- pyFunction("sum")
#' pySum(1:3)
#' }
#  ---------------------------------------------------------
pyFunction <- function(key){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    cfun <- sprintf(callFun, key)
    fun <- eval(parse(text=cfun))
    class(fun) <- "pyFunction"
    attr(fun, "name") <- key
    fun
}

print.pyFunction <- function(x, ...) pyExecp(attr(x, "name"))
