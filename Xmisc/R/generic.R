
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Sat Jan 11 22:01:40 EST 2014 -0500 (Week 01)
## 
## 
## Reference: 
## 
## 
## ************************************************************************


## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------

##setGeneric("popone", function(x,...) standardGeneric("popone"))


##' @title popone-methods
##' @name popone-methods
##' @aliases popone
##' @param obj object
##' @param ... further arguments
##' @param envir the \code{environment} to use
##' @export
##' @docType methods
##' @rdname popone-methods
setGeneric("popone", function(obj,...,envir) {
  if(missing(envir)){
    envir <- sys.frame(sys.parent(2))
  }
  res <- standardGeneric("popone")
  ## logme(deparse(substitute(obj)))
  ## logme(res,"popone")
  assign(deparse(substitute(obj)),res$list,envir=envir)
  return(res$one)
})


##' Remove and return an element from an R object, given name or index.
##' Default: the last element
##' Note: the R object is altered in place, like Python List pop() 
##' 
##' @title Remove and return an element from a list
##' @param obj an R object
##' @param x name or index
##' @param warn logical, whether warn at error
##' @param error logical, whether stop at error
##' @return the element and the R object is altered in place
##' @author Xiaobei Zhao
##' @examples
##' ll <- list(1,2,3,a=4,b=5,6,7,c=8,9)
##' popone(ll,"a")
##' popone(ll,1) ## remove the 1st in ll
##' popone(ll,1) ## remove the next (2nd of the origin)
##' try(popone(list(),"a"))
##' 
setMethod("popone", c("list"), function(
  obj, x, warn=TRUE, error=TRUE){
  'Remove the one at the given index/position (or with the given name) in the list, and return it. If no index is specified, popone(obj) removes and returns the last one in the list.'
  ret <- list(one=NULL,list=obj)
  if (!length(obj)){
    msg <- "list must not be empty to pop."
    if (error) {stop(msg)} else {if (warn) message(msg)}
  }
  if (missing(x)){
    x <- length(obj)
  } 
  .x <- x
  if (is.character(x)){
    x <- which(names(obj) %in% x)
    if (length(x)){
      x <- max(x) # pop the last
    } else {
      msg <- sprintf('popone | no such one: "%s"',.x)
      if (error) {stop(msg)} else {if (warn) message(msg)}
    }
  }
  if (length(x)!=1){
    msg <- "popone | x must be of length 1"
    if (error) {stop(msg)} else {if (warn) message(msg)}
  } else {
    ret <- list(one=obj[[x]],list=obj[-x])
    ## obj <<- obj[-x] #assign
    ## .name <- deparse(substitute(obj))
    ## .envir <- sys.frame(sys.parent(2))
    ## logme(substitute(obj)) # the argname in setGeneric
    ## logme(.name)
    ## logme(.envir)
    ## assign(.name,obj[-x],envir=.envir)
  }
  return(ret)
})


## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------

##' @title removeone-methods
##' @name removeone-methods
##' @aliases removeone
##' @param obj object
##' @param ... further arguments
##' @param envir the \code{environment} to use
##' @export 
##' @docType methods
##' @rdname removeone-methods
setGeneric("removeone",  function(obj,...,envir=sys.frame(sys.parent(2))) {
  ## .envir <- sys.frame(sys.parent(2))
  ## logme(.envir)
  if(missing(envir)){
    envir <- sys.frame(sys.parent(2))
  }
  res <- standardGeneric("removeone")
  assign(deparse(substitute(obj)),res,envir=envir)
  return()
})


##' Remove an element from an R object
##' Note: the R object is altered in place 
##' 
##' @title Remove an element from an R object
##' @param obj an R object
##' @param x an element
##' @param warn logical, whether warn at error
##' @param error logical, whether stop at error
##' @return NULL and the R object is altered in place
##' @author Xiaobei Zhao
##' @examples
##' ll=list(1,2,3,a=4,b=5,6,7,c=8,9)
##' removeone(ll,3)
##' 
setMethod("removeone", c("list"), function(
  obj, x, warn=TRUE, error=TRUE){
  'Remove the first one from the list whose element is x. It is an error if there is no such one.'
  if (!length(obj)){
    msg <- 'list must not be empty to `remove`.'
    if (error) {stop(msg)} else {if (warn) message(msg)}
  }
  which.x <- which(obj %in% x)
  if (!length(which.x)){
    msg <- sprintf("removeone | Not found: %s",x)
    if (error) {stop(msg)} else {if (warn) message(msg)}
  }
  which.x <- min(which.x) # the last
  ## obj <<- obj[-which.x]
  ret <- obj[-which.x]
  return(ret)
})


## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------

##' @title popmany-methods
##' @name popmany-methods
##' @aliases popmany
##' @param obj object
##' @param ... further arguments
##' @param envir the \code{environment} to use
##' @export
##' @docType methods
##' @rdname popmany-methods
setGeneric("popmany", function(obj,...,envir) {
  if(missing(envir)){
    envir <- sys.frame(sys.parent(2))
  }
  res <- standardGeneric("popmany")
  assign(deparse(substitute(obj)),res$list,envir=envir)
  return(res$many)
})


##' Remove and return many elements from a list
##' 
##' @title Remove and return many elements from a list
##' @param obj list
##' @param x indexes
##' @return the elements and the R object is altered in place
##' @author Xiaobei Zhao
setMethod("popmany", c("list"), function(obj,x){
  ret <- list(many=obj[x],list=obj[-x])
  return(ret)
})



## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------
