
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Thu Jan 09 07:52:27 EST 2014 -0500 (Week 01)
## 
## 
## Reference: 
## 
## 
## ************************************************************************



## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------

`.` <- function (..., .env=parent.frame())
{
  structure(as.list(match.call()[-1]), env=.env, class="quoted")
}


`..` <- function (..., .env=sys.parent(2))
{
  get(deparse(substitute(...)), envir=.env)
}


## ------------------------------------------------------------------------
## logger
## ------------------------------------------------------------------------

##' Define log levels
##'
##' 
##' @title Define log levels
##' @return the defined loglevels
##' @author Xiaobei Zhao
get_loglevel <- function(){
  c('NULL','INFO', 'DEBUG', 'WARNING', 'ERROR', 'CRITICAL')
}

##' General test of an object being interpretable as \code{loglevel}
##'
##' 
##' @title General test of an object being interpretable as loglevel
##' @param x object
##' @param loglevels the defined loglevels
##' @return logical
##' @author Xiaobei Zhao
is.loglevel <- function(
  x,
  loglevels=get_loglevel()
  )
{
  x %in% loglevels
}

##' Coerces an object to \code{loglevel}
##'
##' 
##' @title Coerces an object to loglevel
##' @param x object
##' @param loglevels the defined loglevels
##' @return loglevel
##' @author Xiaobei Zhao
as.loglevel <- function(
  x,
  loglevels=get_loglevel()
  )
{
  if (is.null(x)){
    x <- 'NULL'
  }
  if (!is.loglevel(x,loglevels=loglevels)){
    stop('Invalid loglevel')
  }
  if ( is.factor(x) ){
    if (identical(levels(x),loglevels)){
      return(x)
    }
  }
  x <- factor(x,levels=loglevels)
  
  return(x)
}
  



## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------

##' General test of a class name being activeBindingFunction
##'
##' 
##' @title General test of a class name being activeBindingFunction
##' @param cls the name of a class
##' @return logical
##' @author Xiaobei Zhao
is.activeBindingFunction <- function(cls){
  "activeBindingFunction" == cls
}

##' General test of a class being uninitializedField
##'
##' 
##' @title General test of a class being uninitializedField
##' @param x class
##' @return logical
##' @author Xiaobei Zhao
##' @seealso \code{methods::ReferenceClasses}
is.uninitializedField <- function(x){
  "uninitializedField" == class(x)
}



## ------------------------------------------------------------------------
## time
## ------------------------------------------------------------------------

##' Get system epoch
##'
##' 
##' @title Get system epoch
##' @return numeric
##' @author Xiaobei Zhao
Sys.Epoch <-
  function()
{
  .digits <- options()$digits
  options(digits=22)
  ret <- as.numeric(as.POSIXct(Sys.time()))
  options(digits=.digits)
  return(ret)
}


## ------------------------------------------------------------------------
## generic
## ------------------------------------------------------------------------


##' @title getone-methods
##' @name getone-methods
##' @aliases getone
##' @param obj object
##' @param ... further arguments
##' @export
##' @docType methods
##' @rdname getone-methods
setGeneric("getone", function(obj, ...) standardGeneric("getone"))


##' Get an element by index or name from a list
##'
##' 
##' @title Get an element by index or name from a list
##' @param obj list
##' @param x index or name of an element
##' @param safe whether allow to get an element by index,
##' in case this element is after any named element
##' @param msg whether print a message
##' @return list
##' @author Xiaobei Zhao
##' @examples
##' ll <- list(11,22,33,a=44,b=55,66,77,c=88,99)
##' getone(ll,numeric())
##' getone(ll,"a")
##' getone(ll,"c")
##' getone(ll,1)
##' getone(ll,7)
##' getone(ll,7,safe=FALSE)
##' 
setMethod("getone", c("list"), function(
  obj,x,safe=TRUE,msg=TRUE){
  ## logme(obj,'setMethod | getone','DEBUG')
  ## logme(x,'setMethod | getone','DEBUG')
  
  .idx <- valid.arg.index(obj,x,safe=safe)
  if (length(.idx)){
    ret <- obj[[.idx]]
  } else {
    ret <- NULL
  }
  if (safe){
    .idx0 <- valid.arg.index(obj,x,safe=FALSE)
    if (!length(.idx) & length(.idx0)){
      .names <- names(obj)
      if (msg) {
        message(sprintf('List | getone | safe: 
idx=%s; .names=c(%s)',.idx0, paste('"',.names,'"',sep='',collapse=',')))
      }
    }
  }
  
  ## logme(ret,'setMethod | getone','DEBUG')
  return(ret)
})



## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------

##' Check validity of an index of a \code{list} object. 
##'
##' @title Check validity of an index
##' @param obj list
##' @param x index or name of an element
##' @param safe whether safe if an index is higher than the one of any named element
##' @author Xiaobei Zhao
##' @return numeric, the index. 
##' Return numeric(0) if the name or the index does not exist or when the
##' index is invalid. 
##' If \code{safe} is FALSE, any index is valid; if \code{safe} is TRUE,
##' an index is invalid when the indexed element is positionally after
##' another named element.  
##' 
##' @examples
##' ll <- list(11,12,13,a=14,b=15,16,17,c=18,19)
##' valid.arg.index(ll,-1) # non-existing index
##' valid.arg.index(ll,0)  # non-existing index
##' valid.arg.index(ll,1)  # valid index
##' valid.arg.index(ll,2)  # valid index
##' valid.arg.index(ll,5)  # invalid index
##' valid.arg.index(ll,10) # non-existing index
##' valid.arg.index(ll,"a")# valid name
##' valid.arg.index(ll,"e")# non-existing name
##' valid.arg.index(ll,5,safe=FALSE) # still return the index
##' 
valid.arg.index <-
  function(obj,x,safe=TRUE)
{
  .names <- names(obj)
  if (is.character(x)){
    if (x %in% .names) {
      x <- which(.names %in% x)[1] # only the first match
      ret <- TRUE
    } else {
      ret <- FALSE
    }
  } else if (is.numeric(x) | is.integer(x)){
    if (!length(x)){
      ret <- FALSE
    } else{
      ret <- (x>=1 & x<=length(obj))
      if (safe & length(.names)){
        max.idx <- min(which(.names!=""))-1
        ret <- (x>=1 & x<=max.idx)
      }
    }
  } else {
    stop("valid.arg.index | x must be character/numeric/integer")
  }
  if (ret) {
    ret <- x
  } else {
    ret <- numeric()
  }
  return(ret)
}


