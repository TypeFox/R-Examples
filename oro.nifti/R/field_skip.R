#' @name field_skip-methods
#' @title Extract Image Attribute \code{field_skip}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{field_skip} field.  
#' @description Methods that act on the \code{field_skip} field in the
#' NIfTI/ANALYZE header.
#' @rdname field_skip-methods
#' @aliases field_skip-methods, field_skip
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
setGeneric("field_skip", function(object) standardGeneric("field_skip"))
#' @rdname field_skip-methods
#' @aliases field_skip,anlz-method
#' @export
setMethod("field_skip", "anlz", function(object) { object@"field_skip" })
#' @rdname field_skip-methods
#' @aliases field_skip<- 
#' @export
setGeneric("field_skip<-", function(object, value) { standardGeneric("field_skip<-") })
#' @rdname field_skip-methods
#' @aliases field_skip<-,anlz-method
#' @export
setMethod("field_skip<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "field_skip" %in% slotNames(object) ){
              object@"field_skip" <- value
            } else {
              warning("field_skip is not in slotNames of object")
            }
            return(object)
          })
#' @rdname field_skip-methods
#' @aliases field.skip,nifti-method
#' @export
setGeneric("field.skip", function(object) standardGeneric("field.skip"))
#' @rdname field_skip-methods
#' @aliases field.skip,anlz-method
#' @export
setMethod("field.skip", "anlz", function(object) { object@"field_skip" })
#' @rdname field_skip-methods
#' @aliases field.skip<- 
#' @export
setGeneric("field.skip<-", function(object, value) { standardGeneric("field.skip<-") })
#' @rdname field_skip-methods
#' @aliases field.skip<-,anlz-method
#' @export
setMethod("field.skip<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "field_skip" %in% slotNames(object) ){
              object@"field_skip" <- value
            } else {
              warning("field_skip is not in slotNames of object")
            }
            return(object)
          })
