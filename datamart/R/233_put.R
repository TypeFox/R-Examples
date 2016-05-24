#' Put a target
#'
#' This generic creates the target at the given location. 
#'
#' @param target    an object of class Target or derived
#' @param where     an object of class Location or derived
#' @param verbose   print diagnostic messages
#' @param ...       additional parameters
#'
#' @export
#' @docType methods
#' @rdname put-methods
setGeneric(
  name="put",
  def=function(target, where, ...){standardGeneric("put")}
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,Target,character-method
setMethod(
  f="put",
  signature=c(target="Target", where="character"),
  definition=function(target, where, ...) {
    if(where==":memory:") 
      newloc <- new("MemoryLocation")
    else if(isTRUE(file.info(where)$isdir)) 
      newloc <- dirloc(path=where)
    else
      stop("cannot interpret 'where' argument '", where, "', please pass ':memory:', a path to an existing directory or an Location object.")
    put(target, newloc, ...)
  }
)
