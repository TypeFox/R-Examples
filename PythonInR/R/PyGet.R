## ----------------------------------------------------------------------------- 
##
##   PyGet always creates a new Python Object if it can't be
##         to R. The new Python Object is stored in the dict __R__.
##
##   PyGet0 creates 0 new Python objects.
##
## -----------------------------------------------------------------------------

#  -----------------------------------------------------------
#  pyGet0
#  ======
#' @title Creates an R representation of an Python object
#'
#' @description The function pyGet0 gets Python objects by name.
#' @param key a string specifying the name of a Python object.
#' @details Primitive data types like bool, int, long, float, str, 
#'          bytes and unicode are returned as R objects. Python tuples, lists,
#'          dictionaries and other Python objects are returned as virtual
#'          Python objects.
#' @note pyGet0 never creates a new Python object.
#' @return Returns the specified Python object converted into an R object if
#'         possible, else a virtual Python object.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' if ( pyIsConnected() ){
#' pyExec("import os")
#' os <- pyGet0("os")
#' os$getcwd()
#' os$sep
#' os$sep <- "Hello Python!"
#' pyExecp("os.sep")
#' }
# -----------------------------------------------------------
pyGet0 <- function(key){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    check_string(key)

    pyClass <- pyType(key)
    
    if (is.null(pyClass)){
        stop(sprintf('"%s" does not exist', key))
    }

    if ( pyClass %in% c("NoneType", "bool", "int", "long", "float", "str", "bytes", "unicode")){
        return(pyGet(key))
    }else if (pyIsCallable(key)){
        return(pyFunction(key))
    }else if ( pyClass == "tuple" ){
        return(pyTuple(key, regFinalizer = FALSE))
    }else if ( pyClass == "list" ){
        return(pyList(key, regFinalizer = FALSE))
    }else if ( pyClass == "dict" ){
        return(pyDict(key, regFinalizer = FALSE))
    }else{
        return(pyObject(key, regFinalizer = FALSE))
    }

}

#  -----------------------------------------------------------
#  pyGet
#  =====
#' @title Gets Python objects by name and transforms them into R objects
#'
#' @description The function pyGet gets Python objects by name and transforms 
#'              them into R objects. 
#' @param key a string specifying the name of a Python object.
#' @param autoTypecast an optional logical value, default is TRUE, specifying
#'        if the return values should be automatically typecasted if possible.
#' @param simplify an optional logical value, if TRUE R converts Python lists 
#'        into R vectors whenever possible, else it translates Python lists 
#'        always into R lists.
#' @details Since any Python object can be transformed into one of the basic 
#'          data types it is up to the user to do so up front. More information
#'          about the type conversion can be found in the README file or at
#'          \url{http://pythoninr.bitbucket.org/}. \cr
#' @note pyGet always returns a new object, if you want to create a R representation
#'       of an existing Python object use pyGet0 instead.
#' @return Returns the specified Python object converted into an R object if
#'         possible, else a virtual Python object.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' ## get a character of length 1
#' pyGet("__name__")
#' ## get a character of length 1 > 1
#' pyGet("sys.path")
#' ## get a list
#' pyGet("sys.path", simplify = FALSE)
#' ## get a PythonInR_List
#' x <- pyGet("sys.path", autoTypecast = FALSE)
#' x
#' class(x)
#'
#' ## get an object where no specific transformation to R is defined
#' ## this example also shows the differnces between pyGet and pyGet0
#' pyExec("import datetime")
#' ## pyGet creates a new Python variable where the return value of pyGet is
#' ## stored the name of the new reference is stored in x$py.variableName.
#' x <- pyGet("datetime.datetime.now().time()")
#' x
#' class(x)
#' x$py.variableName
#' ## pyGet0 never creates a new Python object, objects which can be transformed 
#' ## to R objects are transformed. For all other objects an PythonInR_Object is created.
#' y <- pyGet0("datetime.datetime.now().time()")
#' y
#' ## An important difference is that the evaluation of x always will return the same
#' ## time, the evaluation of y always will give the new time.
# -----------------------------------------------------------
pyGet <- function(key, autoTypecast=TRUE, simplify=TRUE){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    check_string(key)

    pyClass <- pyExecg(sprintf("try:\n\tx=type(%s).__name__\nexcept:\n\tx=None", key))[['x']]
    class(pyClass) <- pyClass

    if (pyOptions("winPython364")){
        x <- tryCatch(pyGetPoly(key, autoTypecast, simplify, pyClass), 
              error = function(e){
                  err <- e$message
                  class(err) <- "errorMessage"
                  return(err)})
        cat(pyGetSimple("__getStdout()")) ## print stdout
        if (class(x) == "errorMessage"){
            msg <- makeErrorMsg()
            if (!is.null(msg)) stop(msg)
            cat(x)
        }
    }else{
        x <- pyGetPoly(key, autoTypecast, simplify, pyClass)
    }

    return(pyTransformReturn(x))
}

pyGetSimple <- function(key, autoTypecast=TRUE, simplify = TRUE){
  .Call("PythonInR_Run_String", key, 258L, autoTypecast, FALSE, FALSE, 1L, simplify)
}

pyGetPoly <- function(key, autoTypecast, simplify, pyClass){
  pyGetSimple(key, autoTypecast, simplify)
}

setGeneric("pyGetPoly")

setClass("PrVector")
setMethod("pyGetPoly", signature(key = "character", autoTypecast = "logical", simplify = "logical", pyClass = "PrVector"),
          function(key, autoTypecast, simplify, pyClass){
    x <- pyExecg(sprintf("x = %s.toDict()", key), autoTypecast = autoTypecast, simplify = simplify)[['x']]
    y <- setNames(x[['vector']], x[['names']])
    if (is.null(y)) class(y) <- x[['rClass']]
    y
})

setClass("PrMatrix")
setMethod("pyGetPoly", signature(key="character", autoTypecast = "logical", simplify = "logical", pyClass = "PrMatrix"),
          function(key, autoTypecast, simplify, pyClass){
    x <- pyExecg(sprintf("x = %s.toDict()", key), autoTypecast = autoTypecast, simplify = simplify)[['x']]
    M <- do.call(rbind, x[['matrix']])
    if ( is.null(M) ) M <- matrix()
    rownames(M) <- x[['rownames']]
    colnames(M) <- x[['colnames']]
    return(M)
})

setClass("ndarray")
setMethod("pyGetPoly", signature(key="character", autoTypecast = "logical", simplify = "logical", pyClass = "ndarray"),
          function(key, autoTypecast, simplify, pyClass){
    x <- pyExecg(sprintf("x = %s.tolist()", key), autoTypecast = autoTypecast, simplify = simplify)[['x']]
    return( do.call(rbind, x) )
})

setClass("PrDataFrame")
setMethod("pyGetPoly", signature(key="character", autoTypecast = "logical", simplify = "logical", pyClass = "PrDataFrame"),
          function(key, autoTypecast, simplify, pyClass){
    x <- pyExecg(sprintf("x = %s.toDict()", key), autoTypecast = autoTypecast, simplify = simplify)[['x']]
    df <- as.data.frame(unname(x['data.frame']), stringsAsFactors=FALSE)
    rownames(df) <- x[['rownames']]
    return(df)
})

setClass("DataFrame")
setMethod("pyGetPoly", signature(key="character", autoTypecast = "logical", simplify = "logical", pyClass = "DataFrame"),
          function(key, autoTypecast, simplify, pyClass){
    x <- pyExecg(sprintf("x = %s.to_dict()", key), autoTypecast = autoTypecast, simplify = simplify)[["x"]]
    return( as.data.frame(x, stringsAsFactors=FALSE) )
})

