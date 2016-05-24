# ~ try() but catches warnings : 
##' Catch *and* save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
tryCatch.W.E <- function(expr) {
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
   				     warning = w.handler),
       warning = W)
}

## meeting the width.cutoff in deparse :-(
DEPARSE <- function(expr) {
  paste(deparse(expr),collapse="")
}

CBIND <- function(x,y) UseMethod("CBIND") 
CBIND.matrix <- function(x,y) cbind(x,as.matrix(y))
CBIND.Matrix <- function(x,y) cBind(x,y)

RBIND <- function(x,y) UseMethod("RBIND") 
RBIND.matrix <- function(x,y) rbind(x,as.matrix(y))
RBIND.Matrix <- function(x,y) rBind(x,y)

overcat <- function(msg, prevmsglength) {
  msglength <- nchar(msg)
  if (prevmsglength>0) {cat("\r")}    
  cat(msg)
  return(msglength)
}

