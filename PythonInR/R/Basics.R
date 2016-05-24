# ------------------------------------------------------------------------------ 
#
#   Basics
#
# ------------------------------------------------------------------------------

#  -----------------------------------------------------------
#  pyDir
#  =====
#' @title Convenience function to call the Python function \strong{dir}
#'
#' @description A convenience function to call the Python function \strong{dir}.
#' @param objName an optional string specifying the name of the Python object.
#' @return Returns the list of names in the global scope, if no object name is 
#'         provided, otherwise a list of valid attributes for the specified object.
#' @details The Python function dir is similar to the R function ls.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' pyDir()
#' pyDir("sys")
#  -----------------------------------------------------------
pyDir <- function(objName=NULL){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    if ( is.null(objName) ){
        cmd <- '__tmp__=dir()'
    }else{
        check_string(objName)
        cmd <- sprintf('__tmp__=dir(%s)', objName)
    }
    # Looks redundant but is necessary because of the different namespaces!
    # If I would call pyExecg("x=dir()"), I would only get the elements of
    # the temporary namespace.
    pyExec(cmd)
    retv = pyExecg("__tmp__=__tmp__")[[1]]
    pyExec("__tmp__=None;del(__tmp__)")
    retv
}

#  ----------------------------------------------------------- 
#  pyHelp
#  ======
#' @title Convenience function to access the Python \strong{help} system
#'
#' @description a convenience function to access the Python \strong{help} system.
#' @param topic a string specifying name or topic for which help is sought.
#' @return Prints the help to the given string.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' pyHelp("abs")
#  -----------------------------------------------------------
pyHelp <- function(topic){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    check_string(topic)
    pyExecp(sprintf("help('%s')", topic))
}

#  -----------------------------------------------------------
#  pyType
#  ======
#' @title Convenience function to call the Python function \strong{type}
#'
#' @description Convenience function to call the Python function \strong{type}.
#' @param objName a string specifying the name of the Python object.
#' @return The type of the specified object as character on success, NULL otherwise.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' pyExec("x = dict()")
#' pyType("x")
#  -----------------------------------------------------------
pyType <- function(objName){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    check_string(objName)
    cmd <- 'try:\n\tx = type(%s).__name__\nexcept:\n\tx = None'
    retVal <- tryCatch({pyExecg(sprintf(cmd, objName))[['x']]},
        error=function(e){
            errMsg <- sprintf('pyType("%s")\n  >>> type(%s)\n  SyntaxError: invalid syntax',
                                objName, objName)
            class(errMsg) <- "errorMessage"
            return(errMsg)
        })
    if (class(retVal) == "errorMessage") stop(retVal, call.=FALSE)
    return(retVal)
}
