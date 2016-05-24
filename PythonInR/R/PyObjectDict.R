##  
##  pyDict
##
##
##  PythonInR_Dict
##      class definition

PythonInR_Dict <-
  R6Class(
    "PythonInR_Dict",
    portable = TRUE,
    inherit = PythonInR_Object,
    public = list(
      print = function() pyExecp(self$py.variableName),
      clear = function(){
          cable <- sprintf("%s.clear", self$py.variableName)
          pyCall(cable)},
      copy = function(){
          cable <- sprintf("%s.copy", self$py.variableName)
          pyCall(cable)},
      fromkeys = function(seq, value){
          cable <- sprintf("%s.fromkeys", self$py.variableName)
          if (missing(value)) pyCall(cable, list(seq))
          else pyCall(cable, list(seq, value))},
      get = function(key, default){
          cable <- sprintf("%s.get", self$py.variableName)
          if (missing(value)) pyCall(cable, list(key))
          else pyCall(cable, list(key, default))},
      has_key = function(key){
          cable <- sprintf("%s.has_key", self$py.variableName)
          pyCall(cable, list(key))},
      items = function(){
          cable <- sprintf("%s.items", self$py.variableName)
          pyCall(cable)},
      keys = function(key){
          cable <- sprintf("%s.keys", self$py.variableName)
          pyCall(cable)},
      pop = function(key, default){
          cable <- sprintf("%s.pop", self$py.variableName)
          if (missing(value)) pyCall(cable, list(key))
          else pyCall(cable, list(key, default))},
      popitem = function(){
          cable <- sprintf("%s.popitem", self$py.variableName)
          pyCall(cable)},
      setdefault = function(key, default){
          cable <- sprintf("%s.setdefault", self$py.variableName)
          if (missing(value)) pyCall(cable, list(key))
          else pyCall(cable, list(key, default))},
      update = function(dict){
          cable <- sprintf("%s.update", self$py.variableName)
          pyCall(cable, list(dict))},
      values = function(){
          cable <- sprintf("%s.values", self$py.variableName)
          pyCall(cable)},
      viewitems = function(){
          cable <- sprintf("%s.viewitems", self$py.variableName)
          pyCall(cable)},
      viewkeys = function(){
          cable <- sprintf("%s.viewkeys", self$py.variableName)
          pyCall(cable)},
      viewvalues = function(){
          cable <- sprintf("%s.viewvalues", self$py.variableName)
          pyCall(cable)}
    ))

PythonInR_DictNoFinalizer <-
    R6Class("PythonInR_Dict",
            portable = TRUE,
            inherit = PythonInR_Dict,
            public = list(
                initialize = function(variableName, objectName, type) {
                    if (!missing(variableName)) self$py.variableName <- variableName
                    if (!missing(objectName)) self$py.objectName <- objectName
                    if (!missing(type)) self$py.type <- type
                }
            ))

`[.PythonInR_Dict` <- function(x, i){
    slice <- deparse(i)
    pyGet(sprintf("%s[%s]", x$py.variableName, slice))
}

`[<-.PythonInR_Dict` <- function(x, i, value){
    if (length(i) > 1) class(i) <- "tuple"
    success <- .Call("r_set_py_dict", x$py.variableName, i, value)
    x
}

#  ---------------------------------------------------------
#  pyDict
#  ======
#' @title Create a virtual Python dictionary
#'
#' @description The function pyDict creates a virtual Python object 
#'              of type PythonInR_Dict.
#' @param key a character string giving the name of the Python object.
#' @param value if a value is provided, a new Python dictionary is created based 
#'              on the value. Therefore allowed values of value are named lists and 
#'              names vectors.
#' @param regFinalizer a logical indicating if a finalizer should be
#'                     be registered, the default value is TRUE.
#' @details If no value is provided a virtual Python dict for an existing
#'          Python object is created. If the value is NULL, an empty 
#'          virtual Python object for an empty dict is created.
#'          If the value is a named vector or named list, a new Python
#'          object based on the vector or list is created.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' if ( pyIsConnected() ){
#' pyExec('myPyDict = {"a":1, "b":2, "c":3}')
#' ## create a virtual Python dictionary for an existing dictionary
#' myDict <- pyDict("myPyDict")
#' myDict["a"]
#' myDict["a"] <- "set the key"
#' myDict
#' ## allowed keys are
#' myDict['string'] <- 1
#' myDict[3L] <- "long"
#' myDict[5] <- "float"
#' myDict[c("t", "u", "p", "l", "e")] <- "tuple"
#' myDict
#' ## NOTE: Python does not make a difference between a float key 3 and a long key 3L
#' myDict[3] <- "float"
#' myDict
#' ## create a new Python dict and virtual dict
#' myNewDict <- pyDict('myNewDict', list(p=2, y=9, r=1))
#' myNewDict
#' }
#  ---------------------------------------------------------
pyDict <- function(key, value, regFinalizer = TRUE){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    check_string(key)

    if (!missing(value)){
        ## create a new object
        if (is.null(value)){
          pyExec(sprintf("%s = dict()", key))
        }else{
          pySetSimple(key, value)
        }
    }

    if (!pyVariableExists(key))
        stop(sprintf("'%s' does not exist in the global namespace",
             key))
    vIsDict <- pyGet(sprintf("isinstance(%s, dict)", key))
    if (!vIsDict)
        stop(sprintf("'%s' is not an instance of dict", key))

    if (regFinalizer){
      py_dict <- PythonInR_Dict$new(key, NULL, "dict")
    }else{
      py_dict <- PythonInR_DictNoFinalizer$new(key, NULL, "dict")
    }
    return(py_dict)
}
