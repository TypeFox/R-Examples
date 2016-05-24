# ------------------------------------------------------------------------------ 
#
#   PyCall
#
# ------------------------------------------------------------------------------

#  -----------------------------------------------------------
#  pyCall
#  ======
#' @title Call a callable Python object from within R
#'
#' @description Call a callable Python object from within R.
#' @param callableObj a character string giving the name of the desired callable 
#'                    Python object.
#' @param args an optional list of arguments passed to the callable. 
#' @param kwargs an optional list of named arguments passed to the callable.
#' @param autoTypecast an optional logical value, default is TRUE, specifying
#'        if the return values should be automatically typecasted if possible.
#' @param simplify an optional logical value, if TRUE, R converts Python lists 
#'                 into R vectors whenever possible, else it translates Python 
#'                 lists always to R lists.
#' @return Returns the result of the function call, converted into an R object.
#' @details The args and kwargs are transformed to Python variables by the 
#'          default conversion. More information about the type conversion can 
#'          be found in the vignette.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' pyCall("sum", args=list(1:3))
#' 
#' ## define a new function with the name fun
#' pyExec('
#' def fun(**kwargs):
#'     return([(key, value) for key, value in kwargs.items()])
#' ')
#' pyCall("fun", kwargs=list(a=1, f=2, x=4))
#  -----------------------------------------------------------
pyCall <- function(callableObj, args=NULL, kwargs=NULL, autoTypecast=TRUE, simplify=TRUE){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    check_string(callableObj)
    
    if (pyOptions("winPython364")){
        returnValue <- try(.Call("py_call_obj", callableObj, args, kwargs, simplify, autoTypecast), 
                           silent=TRUE)
        msg <- makeErrorMsg()
        if (!is.null(msg)) stop(msg)
    }else{
        returnValue <- .Call("py_call_obj", callableObj, args, kwargs, simplify, autoTypecast)
    }
    return(pyTransformReturn(returnValue))
}
