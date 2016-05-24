#' Xdata -- A class representing a data source
#' 
#' Most methods of the class are abstract,
#' however the \code{show}, \code{print},
#' \code{queries} methods will usually not
#' need to be redefined.
#'
#' The \code{query} method is defined for
#' character resource arguments. It is tried
#' to transform the argument to an object; if
#' that succeeds, \code{query} is called again.
#' Derived methods that are also interpreting
#' resource as character should first call
#' this method via \code{callNextMethod}.
#'
#' @name Xdata-class
#' @rdname Xdata-class
#' @exportClass Xdata
setClass(Class="Xdata", representation=representation())
setClass(Class="EmptySet", representation=representation(), contains="Xdata")

#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,Xdata,character-method
setMethod(
  f="query",
  signature=c(self="Xdata", resource="character"),
  definition=function(self, resource, ...) {
     stop(sprintf("Invalid resource '%s' specified for data object '%s'", resource, class(self)))
  }
)

#' @rdname queries-methods
#' @name queries
#' @docType methods
#' @export 
#' @aliases queries queries,Xdata-method
setMethod(
  f="queries",
  signature="Xdata",
  definition=function(self) {
    proc_one <- function(md) tryCatch(attr(md, "target")[["resource"]], error=function(e) NULL) # md is of class MethodDefinition
    md_list <- findMethods("query", classes=c(class(self), names(getClass(class(self))@contains)))
    res <- sapply(md_list, proc_one)
    res <- res[res!="character"]
    names(res) <- NULL
    return(res)
  }
)

#' Show Method for Xdata classes
#'
#' The \code{show} method for the Xdata class has been adapted to display the class name.
#' Some inherited classes such DirectoryLocation or Blogger override this default definition.
#'
#' @param object    Xdata object
#' @rdname show-methods
#' @name show
#' @export
#' @docType methods
#' @aliases show show,Xdata-method
setMethod(
  f="show",
  signature="Xdata",
  definition=function(object) cat(sprintf("<object of class %s>\n", class(object)))
)

#' @rdname meta-methods
#' @name meta
#' @export
#' @docType methods
#' @aliases meta meta,Xdata-method
setMethod(
  f="meta",
  signature=c("Xdata"),
  definition=function(self, ...) data.frame() 
)

#' @rdname dependencies-methods
#' @name dependencies
#' @export
#' @docType methods
#' @aliases dependencies dependencies,Xdata-method
setMethod(
  f="dependencies",
  signature=c("Xdata"),
  definition=function(self) return(list())
)
