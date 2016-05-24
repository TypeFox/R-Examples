#' @name start_field-methods
#' @title Extract Image Attribute \code{start_field}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{start_field} field.  
#' @description Methods that act on the \code{start_field} field in the
#' NIfTI/ANALYZE header.
#' @rdname start_field-methods
#' @aliases start_field-methods, start_field
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
setGeneric("start_field", function(object) standardGeneric("start_field"))
#' @rdname start_field-methods
#' @aliases start_field,anlz-method
#' @export
setMethod("start_field", "anlz", function(object) { object@"start_field" })
#' @rdname start_field-methods
#' @aliases start_field<- 
#' @export
setGeneric("start_field<-", function(object, value) { standardGeneric("start_field<-") })
#' @rdname start_field-methods
#' @aliases start_field<-,anlz-method
#' @export
setMethod("start_field<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "start_field" %in% slotNames(object) ){
              object@"start_field" <- value
            } else {
              warning("start_field is not in slotNames of object")
            }
            return(object)
          })
