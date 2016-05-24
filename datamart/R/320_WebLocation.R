#' Web location
#'
#' Read-only web folder
#'
#'
#' The queries and meta methods do not work for this class.
#' 
#' @examples
#' getSlots("WebLocation")
#'
#' @seealso \code{\link{webloc}}
#'
#' @name WebLocation-class
#' @rdname WebLocation-class
#' @exportClass WebLocation
setClass(
    Class="WebLocation", 
    representation=representation(path="character"), 
    contains="Location",
    validity=function(object) if(!isTRUE(RCurl::url.exists(object@path))) stop("invalid path argument, seems not to be a valid URL")
)

#' Constructor for a WebLocation object
#'
#' This function returns a WebLocation object that can be used as a read-only folder.
#' 
#' @param path    character, pointing to an existing web directory. Required.
#' @param clss    character, optional class name. Default is "WebLocation".
#'
#' @rdname WebLocation-class
#' @export
webloc <- function(path, clss="WebLocation") new(clss, path=path)

#' @rdname show-methods
#' @name show
#' @export
#' @docType methods
#' @aliases show show,WebLocation-method
setMethod(
  f="show",
  signature="WebLocation",
  definition=function(object) cat(sprintf("<Remote Directory @ %s>\n", object@path))
)

#' @rdname as.character-methods
#' @name as.character
#' @export
#' @docType methods
#' @aliases as.character,WebLocation-method
setMethod(
  f="as.character",
  signature="WebLocation",
  definition=function(x) x@path
)

#' @rdname meta-methods
#' @name meta
#' @export
#' @docType methods
#' @aliases meta meta,WebLocation-method
setMethod(
  f="meta",
  signature="WebLocation",
  definition=function(self) NA
)

#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,WebLocation,character-method
setMethod(
    f="query",
    signature=c(self="WebLocation", resource="character"),
    definition=function(self, resource, verbose=getOption("verbose"), extract.fct=readLines, ...) {
        if(verbose) cat("constructing URL..\n")
        uri <- file.path(self@path, resource)
        if(!isTRUE(RCurl::url.exists(uri))) stop("invalid web location: ", uri)
        if(verbose) cat("downloading URL..\n")
        extract.fct(uri)
    }
)
