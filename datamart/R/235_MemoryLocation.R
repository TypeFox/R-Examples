#' S4 base class to represent output in Memory
#'
#' The \code{MemoryLocation2} class represents the place in the memory of
#' the current R process.
#' 
#' @examples
#' getSlots("MemoryLocation2")
#'
#' @seealso \code{\link{memloc}}
#'
#' @name MemoryLocation2-class
#' @rdname MemoryLocation2-class
#' @exportClass MemoryLocation2
setClass(Class="MemoryLocation2", representation=representation(env="environment"), contains="Location")

#' Constructor for MemoryLocation2 objects
#'
#' @param ...             initial objects that are contained in the location
#' @param clss            class to construct. Defaults to MemoryLocation2.
#'
#' @export
#' @rdname MemoryLocation2-class
memloc <- function(clss="MemoryLocation2", ...) {
    env <- new.env()
    content <- list(...)
    for(nm in names(content)) assign(nm, content[[nm]], envir=env)
    res <- new(clss, env=env)
}

#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,MemoryLocation2,character-method
setMethod(
  f="query",
  signature=c(self="MemoryLocation2", resource="character"),
  definition=function(self, resource, ...) {
    res <- try(get(resource, envir=self@env), silent=TRUE)
    if(!inherits(res, "try-error")) return(res) else callNextMethod(self=self, resource=resource, ...)
  }
)

#' @rdname queries-methods
#' @name queries
#' @export
#' @docType methods
#' @aliases queries queries,MemoryLocation2-method
setMethod(
  f="queries",
  signature="MemoryLocation2",
  definition=function(self) ls(envir=self@env)
)

#' @rdname meta-methods
#' @name meta
#' @export
#' @docType methods
#' @aliases meta meta,MemoryLocation2-method
setMethod(
  f="meta",
  signature="MemoryLocation2",
  definition=function(self, ...) mem.info(envir=self@env)
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,Target,MemoryLocation2-method
setMethod(
    f="put",
    signature=c(target="Target", where="MemoryLocation2"),
    definition=function(target, where, ...) assign(target@name, target, envir=where@env)
)
