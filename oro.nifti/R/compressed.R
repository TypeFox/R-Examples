#' @name compressed-methods
#' @title Extract Image Attribute \code{compressed}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{compressed} field.  
#' @description Methods that act on the \code{compressed} field in the
#' NIfTI/ANALYZE header.
#' @rdname compressed-methods
#' @aliases compressed-methods, compressed
#' @details See documentation on the ANALYZE and/or NIfTI data standards for
#' more details.
#' @author John Muschelli \email{muschellij2@@gmail.com},\cr
#' Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @references
#' ANALYZE 7.5\cr
#' \url{https://rportal.mayo.edu/bir/ANALYZE75.pdf}\cr
#' NIfTI-1\cr
#' \url{http://nifti.nimh.nih.gov/}
#'
#' @export
setGeneric("compressed", function(object) standardGeneric("compressed"))
#' @rdname compressed-methods
#' @aliases compressed,anlz-method
#' @export
setMethod("compressed", "anlz", function(object) { object@"compressed" })
#' @rdname compressed-methods
#' @aliases compressed<- 
#' @export
setGeneric("compressed<-", function(object, value) { standardGeneric("compressed<-") })
#' @rdname compressed-methods
#' @aliases compressed<-,anlz-method
#' @export
setMethod("compressed<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "compressed" %in% slotNames(object) ){
              object@"compressed" <- value
            } else {
              warning("compressed is not in slotNames of object")
            }
            return(object)
          })
