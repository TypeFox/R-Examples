#' Dependencies of a Xdata object
#'
#' The \code{dependencies} method returns a list of character or
#' list elements that define the resources the object depends on.
#'
#' By default, NULL is returned.
#'
#' Inherited classes should override this method
#' if necessary.
#'
#' @param self    an Xdata object
#'
#' @export
#' @rdname dependencies-methods
#' @docType methods
setGeneric(
  name="dependencies",
  def=function(self){standardGeneric("dependencies")}
)
