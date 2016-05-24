#' Verbose list of resources
#'
#' The \code{meta} method returns a data.frame with meta information 
#' entities available at the location. By default,
#' zero rows are returned.
#'
#' Inherited classes should override this method
#' if necessary.
#'
#' @param self    an Location object
#' @param ...     additional parameters
#'
#' @export
#' @rdname meta-methods
#' @docType methods
setGeneric(
  name="meta",
  def=function(self, ...){standardGeneric("meta")}
)
