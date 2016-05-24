#  -----------------------------------------------------------
#  pyAttach
#  ========
#' @title Attach Python objects to R
#' @description A convenience function to attach Python objects to R.
#'
#' @param what a character vector giving the names of the Python objects, 
#'             which should be attached to R.
#' @param env the environment where the virtual Python objects are 
#'            assigned to.
#'
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' if ( pyIsConnected() ){
#' pyExec("import os")
#' 
#' ## attach to global
#' ## ----------------
#' ## attach the function getcwd from the module os to R.
#' pyAttach("os.getcwd", .GlobalEnv)
#' os.getcwd
#' os.getcwd()
#' ## attach the string object os.name to R
#' pyAttach("os.name", .GlobalEnv)
#' pyExecp("os.name")
#' os.name
#' ## Since os.name is attached to the globalenv it can be set without using
#' ## the global assignment operator
#' os.name = "Hello Python from R!"
#' pyExecp("os.name")
#' os.name
#' ## Please note if you don't pyAttach to globalenv you have to use 
#' ## the global assignment operator to set the values of the Python objects
#'
#' ## attach to a new environment
#' ## ---------------------------
#' os <- new.env()
#' attach(os, name="python:os")
#' pyAttach(paste("os", pyDir("os"), sep="."), as.environment("python:os"))
#' os.sep
#' os.sep = "new sep" ## this doesn't changes the value in Python but only 
#'                    ## assigns the new variable os.sep to globalenv
#' os.sep
#' .GlobalEnv$`os.sep`
#' as.environment("python:os")$`os.sep`
#' pyExecp("os.sep")
#' ls()
#' ls("python:os")
#' os.sep <<- "this changes the value in Python"
#' .GlobalEnv$`os.sep`
#' as.environment("python:os")$`os.sep`
#' pyExecp("os.sep")
#' }
# -----------------------------------------------------------
pyAttach <- function(what, env = parent.frame()){
  if ( pyConnectionCheck() ) return(invisible(NULL))
  checkType(environment(), "pyAttach", what="character", env="environment")
  
  for (w in what){
    if (!pyVariableExists(w)) stop(w, " doesn't exist")
  }

  for (o in what){
    po <- o
    spo <- unlist(strsplit(po, split = ".", fixed = TRUE))

    if (length(spo) == 1){
      variableName <- po
      o <- NULL
    }else{
      variableName <- paste(spo[-length(spo)], collapse=".")
      o <- spo[length(spo)] 
    }
    
    if (pyIsCallable(po)){ ## callable functions
      cfun <- sprintf(callFun, po)
      cfun <- eval(parse(text=cfun))
      class(cfun) <- "pyFunction"
      attr(cfun, "name") <- po
      retVal <- assign(po, cfun, envir=env)
    }else{ ## active binding functions
      if (is.null(o)){
        afun <- sprintf(activeFun0, variableName, variableName)
      }else{
        afun <- sprintf(activeFun, variableName, o, o, variableName)
      }
      retVal <- makeActiveBinding(po, eval(parse(text=afun)), env=env)
    }
  }
  return(invisible(env))
}

## new version
## pyPolluteSearchPath <- function(variableName, exports, env=parent.env(environment())){
##    ns <-makeNamespace(sprintf("python:%s", variableName))
##    pydir <- pyDir(variableName)
## 
##    actBind <- NULL
##    for (o in pydir){
##        po <- paste(c(variableName, o), collapse=".")
##        if (pyIsCallable(po)){ # callable functions
##            cfun <- sprintf(callFun, variableName, o, variableName, o, variableName, o)
##            assign(po, eval(parse(text=cfun)), envir=ns)
##        }else{ # active binding functions
##            afun <- sprintf(activeFun, variableName, o, o, variableName)
##            makeActiveBinding(po, eval(parse(text=afun)), env=ns)
##            actBind <- c(po, actBind)
##        }
##    }
## 
##    if (is.null(exports)){
##        namespaceExport(ns, ls(ns))
##    }else{
##        namespaceExport(ns, exports)
##    }
##    env <- attachNamespace(ns)
##    for (acb in actBind) unlockBinding(acb, env)
##    invisible(NULL)
## }

## pyPolluteSearchPath <- function(variableName, exports){
##     ## http://r.789695.n4.nabble.com/Active-bindings-in-attached-environments-td920310.html
##     pydir <- pyDir(variableName)
##     py <- new.env()
##     for (o in pydir){
##         po <- paste(c(variableName, o), collapse=".")
##         if (pyIsCallable(po)){ # callable functions
##             cfun <- sprintf(callFun, variableName, o, variableName, o, variableName, o)
##             assign(po, eval(parse(text=cfun)), env=py)
##         }else{ # active binding functions
##             afun <- sprintf(activeFun, variableName, o, o, variableName)
##             makeActiveBinding(po, eval(parse(text=afun)), env=py)
##             #lockBinding(po, ns)
##         }
##     }

##     ## (TODO (DONE) is commented out): 
##     ##       add error checks
##     py2 <- new.env()
##     attach( py2, pos = 2L , name=sprintf("python:%s", variableName))
##     .Internal(importIntoEnv( as.environment(2L), ls(py), py, ls(py) ))
## }
