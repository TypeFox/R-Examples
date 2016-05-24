##  
##  pyTuple
##

##
##  PythonInR_Tuple
##      class definition

PythonInR_Tuple <-
    R6Class("PythonInR_Tuple",
            portable = TRUE,
            inherit = PythonInR_Object,
            public = list(
                print = function() pyExecp(self$py.variableName),
                index = function(x){
                    cable <- sprintf("%s.index", self$py.variableName)
                    pyCall(cable, args=list(x))},
                count = function(x){
                    cable <- sprintf("%s.count", self$py.variableName)
                    pyCall(cable, args=list(x))}
                ))

PythonInR_TupleNoFinalizer <-
    R6Class("PythonInR_Tuple",
            portable = TRUE,
            inherit = PythonInR_Tuple,
            public = list(
                initialize = function(variableName, objectName, type) {
                    if (!missing(variableName)) self$py.variableName <- variableName
                    if (!missing(objectName)) self$py.objectName <- objectName
                    if (!missing(type)) self$py.type <- type
                }
            ))

    

`[.PythonInR_Tuple` <- function(x, i){
    pyGet(sprintf("%s[%s]", x$py.variableName, deparse(i)))
}

`[<-.PythonInR_Tuple` <- function(x, i, value){
    stop("'tuple' object does not support item assignment", call. = FALSE)
}

#  ---------------------------------------------------------
#  pyTuple
#  =======
#' @title Creates a virtual Python tuple
#'
#' @description The function pyTuple creates a virtual Python object 
#'              of type PythonInR_Tuple. 
#' @param key a character string giving the name of the Python object.
#' @param value an optional value, allowed values are vectors, lists and NULL.
#' @param regFinalizer a logical indicating if a finalizer should be
#'                     be registered, the default value is TRUE.
#' @details If no value is provided a virtual Python tuple for an existing
#'          Python object is created. If the value is NULL an empty 
#'          virtual Python object for an empty tuple is created.
#'          If the value is a vector or tuple a new Python
#'          object based on the vector or list is created.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' if ( pyIsConnected() ){
#' pyExec('myPyTuple = (1, 2, 5, "Hello R!")')
#' # create a virtual Python tuple for an existing tuple
#' myTuple <- pyTuple("myPyTuple")
#' myTuple[0]
#' tryCatch({myTuple[1] <- "should give an error since tuple are not mutable"},
#'          error = function(e) print(e))
#' myTuple
#' # create a new Python tuple and virtual tuple
#' newTuple <- pyTuple('myNewTuple', list(1:3, 'Hello Python'))
#' newTuple[1]
#' }
#  ---------------------------------------------------------
pyTuple <- function(key, value, regFinalizer = FALSE){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    check_string(key)

    if (!missing(value)){
        if ( !is.vector(value) ) stop("'value' has to be a vector or list")
        if (length(value) < 1) value <- as.list(value)
        class(value) <- "tuple"
        pySetSimple(key, value)       
    }
    
    if (!pyVariableExists(key))
        stop(sprintf("'%s' does not exist in the global namespace",
             key))
    vIsTuple <- pyGet(sprintf("isinstance(%s, tuple)", key))
    if (!vIsTuple)
        stop(sprintf("'%s' is not an instance of tuple", key))

    if (regFinalizer){
        py_tuple <- PythonInR_Tuple$new(key, NULL, "tuple")
    }else{
        py_tuple <- PythonInR_TupleNoFinalizer$new(key, NULL, "tuple")
        class(py_tuple) <- class(py_tuple)[-2]
    }
    return(py_tuple)
}
