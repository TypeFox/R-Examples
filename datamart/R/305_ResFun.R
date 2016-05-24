#' ResFunc -- A class representing a "calculated" resource
#' 
#' While the \code{UrlData} and \code{InternalData} classes
#' are examples for accessing actual data, the \code{ResFunc}
#' class represents resources that are not
#' physical data but for example simulated data or
#' data derived from physical data.
#'
#' This class is intended to be used with the
#' \code{Mashup} class.
#'
#' @seealso \code{\link{resfunc}}, \code{\link{resalias}}, \code{\link{datamart}}
#'
#' @name ResFunc-class
#' @rdname ResFunc-class
#' @exportClass ResFunc
setClass(
    Class="ResFunc", 
    representation=representation(
        fun="function", 
        resource="character", 
        arglist="list", 
        depends="list"
    ), 
    contains="Xdata"
)

#' Constructor for ResFunc objects
#'
#' This function creates an ResFunc object. When queried,
#' it returns the result of an function call.
#'
#' @param resource        the name of the resource. Required.
#' @param fun             a function that matches the signature function(self, resource, ...)
#' @param depends         the names of the resources this function depends on.
#' @param clss            name of the class to create. Default ResFunc, must be inherited from this class.
#' @param ...             additional parameters passed to the function when the resource is queried.
#'
#' @export
resfunc <- function(resource, fun, depends=list(), clss="ResFunc", ...) {
    if(!("..." %in% names(formals(fun)))) {
        if(length(setdiff(c("resource", depends), names(formals(fun))))>0) 
            stop(
              "invalid argument list for function resource: need either ellipsis '...' or 'resource' and '", 
              paste(depends, collapse="', '"), 
              "' argument."
            )
    }
    new(clss, resource=resource, fun=fun, arglist=list(...), depends=as.list(depends))
}

#' Phony ResFunc Objects
#'
#' This function creates a ResFunc object that returns an
#' already existing resource. 
#'
#' @param resource        the name of the resource. Required.
#' @param alias_for       name of the resource this object is a proxy for.
#' @param clss            name of the class to create. Default ResFunc, must be inherited from this class.
#'
#' @seealso \code{\link{resfunc}}
#'
#' @export
resalias <- function(resource, alias_for, clss="ResFunc") {
  fun <- function(...) list(...)[[alias_for]]
  resfunc(resource=resource, depends=alias_for, fun=fun, clss=clss)
}

#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,ResFunc,character-method
setMethod(
  f="query",
  signature=c(self="ResFunc", resource="character"),
  definition=function(self, resource, verbose=getOption("verbose"), ...) {
      if(resource==self@resource) {
         arg <- self@arglist
         arg[["verbose"]] <- verbose
         curr_args <- list(...)
         for (nm in names(curr_args)) {
            if(nm!="") arg[[nm]] <- curr_args[[nm]]
            if(is.null(arg[[nm]])) arg[nm] <- list(NULL) # just NULL has no effect
         }
         do.call(self@fun, arg)
      } else {
        if(verbose) cat("trying inherited method..\n")
        callNextMethod(self=self, resource=resource, verbose=verbose, ...)
      }
  }
)

#' @rdname queries-methods
#' @name queries
#' @docType methods
#' @export 
#' @aliases queries queries,ResFunc-method
setMethod(
  f="queries",
  signature="ResFunc",
  definition=function(self) {
    ret <- c(callNextMethod(), self@resource)
    names(ret) <- NULL
    ret <- ret[ret!="character"]
    return(ret)
  }
)


#' @rdname dependencies-methods
#' @name dependencies
#' @export
#' @docType methods
#' @aliases dependencies dependencies,ResFunc-method
setMethod(
  f="dependencies",
  signature=c("ResFunc"),
  definition=function(self) return(self@depends)
)
