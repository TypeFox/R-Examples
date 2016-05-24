#' List resources
#'
#' The \code{queries} method returns a character vector of all
#' defined resources for the given data object. 
#'
#' The default (XData) implementation inspects definitions of the
#' \code{query} method. Inherited classes should override this method
#' if necessary.
#'
#' @param self    an Xdata object
#'
#' @export
#' @docType methods
#' @name queries
#' @rdname queries-methods
setGeneric(
  name="queries",
  def=function(self){standardGeneric("queries")}
)
