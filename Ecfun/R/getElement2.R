getElement2 <- function(object, name=1, default=NA, warn.NULL=TRUE, 
                        envir=list(), returnName){
#       get element of list;  return default if absent
##
## 1.  is.numeric(name)?
##
  if(is.numeric(name)){
    out <- ((name<1) | (length(object)<name))
    if(missing(returnName))
      returnName <- (name == 1)
  } else {
##
## 2.  name not numeric     
##    
    which1 <- which(names(object)==name)
    if(missing(returnName)){
      if(length(which1)==1) {
        returnName <- (which1 == 1)
      } else returnName <- FALSE 
    }
    out <- !(name %in% names(object))
  }
##
## 3.  get object[[name]]
##
  El <- if(out) default else object[[name]]        
##
## 4.  warn.NULL?
##
  if(is.null(El) && warn.NULL && !is.null(default)){
        warning('element ', name, ' is NULL; returning default')
        El <- default
    }
##
## 5.  eval?
##
  if(returnName && is.name(El)){ 
    return(as.character(El))
  } 
  obj.l <- as.list(object)
  Env <- c(obj.l, envir)
  El. <- eval(El, envir=Env) 
  eval(El.)
}

