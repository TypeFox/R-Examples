#' Directory location
#'
#' The show method for the DirectoryLocation class has been adapted to display the path.
#'
#' The as.character method for the DirectoryLocation class returns the path to the directory
#' it represents.
#'
#' The meta method for the DirectoryLocation class returns the output of file.info of the folder.
#' 
#' @examples
#' getSlots("DirectoryLocation")
#'
#' @seealso \code{\link{dirloc}}
#'
#' @name DirectoryLocation-class
#' @rdname DirectoryLocation-class
#' @exportClass DirectoryLocation
setClass(
  Class="DirectoryLocation", 
  representation=representation(path="character"), 
  contains="Location",
  validity=function(object) 
    if(!isTRUE(file.info(object@path)$isdir)) stop("invalid path argument, seems not to be a directory")
)

#' The dirloc function creates a DirectoryLocation object.
#' 
#' @param path    character, pointing to an existing directory. Required.
#' @param clss    character, optional class name. Default is "DirectoryLocation".
#'
#' @rdname DirectoryLocation-class
#' @export
dirloc <- function(path, clss="DirectoryLocation") new(clss, path=path)

#' @rdname show-methods
#' @name show
#' @export
#' @docType methods
#' @aliases show show,DirectoryLocation-method
setMethod(
  f="show",
  signature="DirectoryLocation",
  definition=function(object) cat(sprintf("<Local Directory @ %s>\n", object@path))
)

#' @rdname as.character-methods
#' @name as.character
#' @export
#' @docType methods
#' @aliases as.character,DirectoryLocation-method
setMethod(
  f="as.character",
  signature="DirectoryLocation",
  definition=function(x) x@path
)

#' @rdname meta-methods
#' @name meta
#' @export
#' @docType methods
#' @aliases meta meta,DirectoryLocation-method
setMethod(
  f="meta",
  signature="DirectoryLocation",
  definition=function(self) {
    f <- list.files(path=self@path, include.dirs=FALSE, recursive=TRUE, no..=TRUE)
    fi <- file.info(file.path(self@path, f))
    rownames(fi) <- f
    return(fi)
  }
)

#' @param extract.fct    (DirectoryLocation) which function to use to read file (default \code{readLines})
#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,DirectoryLocation,character-method
setMethod(
    f="query",
    signature=c(self="DirectoryLocation", resource="character"),
    definition=function(self, resource, verbose=getOption("verbose"), extract.fct=readLines, ...) {
        p <- file.path(self@path, resource)
        if(!isTRUE(file.exists(p))) stop("invalid file name: ", p)
        if(verbose) cat("loading file..\n")
        extract.fct(p)
    }
)
