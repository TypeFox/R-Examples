#' Buildable target
#'
#' This is an abstract class for defining buildable targets. 
#'
#' @name Target-class
#' @rdname Target-class
#' @exportClass Target
setClass(Class="Target", representation=representation(name="character"))

#' @rdname show-methods
#' @name show
#' @docType methods
#' @aliases show show,Target-method
#' @export
setMethod(
  f="show",
  signature="Target",
  definition=function(object) cat(sprintf("<object of class %s>\n", class(object)))
)

