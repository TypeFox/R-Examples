#' Request data from data source
#'
#' This generic function is the main interface to
#' the data behind the Xdata layer. The first argument
#' is the data object, the second argument is an identifier
#' (type character), of the resource requested.  
#'
#' Depending on the data object, additional parameter can be
#' provided.
#'
#' @param self        an Xdata object
#' @param resource    an identifier of the resource requested. End-user usually provide character,
#'                    developer use \code{resource} and dispatch on the type.
#' @param verbose     print diagnostic messages, default=FALSE
#' @param ...         additional parameter
#'
#' @export
#' @docType methods
#' @name query
#' @rdname query-methods
setGeneric(
  name="query",
  def=function(self, resource, ...){standardGeneric("query")}
)

