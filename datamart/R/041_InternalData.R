#' A class for querying data()sets
#' 
#' This class allows to query datasets that can be
#' loaded with data(). Only read-only access.
#' 
#' @seealso \code{\link{internalData}}
#'
#' @examples
#' getSlots("InternalData")
#' 
#' @name InternalData-class
#' @rdname InternalData-class
#' @exportClass InternalData
setClass(
    Class="InternalData", 
    representation=representation(name="character", package="character", data_env="environment"), 
    contains="Xdata",
    validity=function(object) if(length(object@data_env)==0) stop("error loading data")
)

#' Constructor for InternalData objects
#'
#' @param name  name of the dataset. Required.
#' @param package name of the package where the dataset is located. Default NULL.
#' @param clss  name of the class to create. Default InternalData, must be inherited from this class.
#'
#' @export
#' @rdname InternalData-class
internalData <- function(name, package, clss="InternalData") {
    e <- new.env()
    do.call("data", list(name=name, package=package, envir=e))
    new(clss, name=name, package=package, data_env=e)
}

#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,InternalData,character-method
# FIXME: constructor's resource is not used
setMethod(
    f="query",
    signature=c(self="InternalData", resource="character"),
    definition=function(self, resource, ...) {
        if(resource %in% ls(envir=self@data_env))
            self@data_env[[resource]]  
        else if(resource==self@resource)
            self@data_env[[self@name]]
        else
            callNextMethod()
    }
)

#' @rdname queries-methods
#' @name queries
#' @export
#' @docType methods
#' @aliases queries queries,InternalData-method
setMethod(
    f="queries",
    signature="InternalData",
    definition=function(self) c(callNextMethod(self=self), ls(envir=self@data_env))
)

#' @rdname meta-methods
#' @name meta
#' @export
#' @docType methods
#' @aliases meta meta,InternalData-method
setMethod(
    f="meta",
    signature="InternalData",
    definition=function(self, ...) mem.info(envir=self@data_env)
)

