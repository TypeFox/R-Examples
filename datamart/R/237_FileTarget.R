#' FileTarget
#' 
#' This class is a decorator for a physical existing file.
#' A common workflow for other targets is to create
#' a temporary file and then call put again with a FileTarget.
#' 
#' @examples
#' getSlots("FileTarget")
#'
#' @seealso \code{\link{filetarget}}
#'
#' @name FileTarget-class
#' @rdname FileTarget-class
#' @exportClass FileTarget
setClass(
    Class="FileTarget", 
    representation=representation(filename="character"),
    contains="Target",
    validity=function(object) 
      if(!file.exists(object@filename)) stop("invalid filename argument, must exist")
)

#' Constructor for FileTarget objects
#'
#' see class FileTarget for details.
#'
#' @param name        name of the Report, default ''
#' @param filename    name of original file
#' @param clss        class name, default 'FileTarget'
#'
#' @return generic
#' @export
#' @rdname FileTarget-class
filetarget <- function(name, filename, clss="FileTarget") 
  new(clss, name=name, filename=filename)


#' @param overwrite    parameter for FileTarget/MdReport -- overwrite existing files? Default TRUE.
#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,FileTarget,DirectoryLocation-method
setMethod(
  f="put",
  signature=c(target="FileTarget", where="DirectoryLocation"),
  definition=function(target, where, overwrite=TRUE) {
    if(target@name=="") stop("Can not put target with empty name into directory location.")
    desti <- file.path(as.character(where), target@name)
    if(file.exists(desti) && overwrite) stop("File '", target@name, "' exists in folder '", as.character(where), "' and overwrite=FALSE")
    file.copy(target@filename, desti)
  }
)

#' @rdname show-methods
#' @name show
#' @export
#' @docType methods
#' @aliases show show,FileTarget-method
setMethod(
  f="show",
  signature="FileTarget",
  definition=function(object) cat(sprintf("<file target %s (%s)>\n", object@name, object@filename))
)


#' @rdname as.character-methods
#' @name as.character
#' @export
#' @docType methods
#' @aliases as.character,FileTarget-method
setMethod(
  f="as.character",
  signature="FileTarget",
  definition=function(x) x@filename
)
