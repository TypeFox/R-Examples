#' @name qoffset_x-methods
#' @title Extract Image Attribute \code{qoffset_x}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{qoffset_x} field.  
#' @description Methods that act on the \code{qoffset_x} field in the
#' NIfTI/ANALYZE header.
#' @rdname qoffset_x-methods
#' @aliases qoffset_x-methods, qoffset_x
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
setGeneric("qoffset_x", function(object) standardGeneric("qoffset_x"))
#' @rdname qoffset_x-methods
#' @aliases qoffset_x,nifti-method
#' @export
setMethod("qoffset_x", "nifti", function(object) { object@"qoffset_x" })
#' @rdname qoffset_x-methods
#' @aliases qoffset_x<- 
#' @export
setGeneric("qoffset_x<-", function(object, value) { standardGeneric("qoffset_x<-") })
#' @rdname qoffset_x-methods
#' @aliases qoffset_x<-,nifti-method
#' @export
setMethod("qoffset_x<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "qoffset_x" %in% slotNames(object) ){
              object@"qoffset_x" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("qoffset_x <-", value))               
            } else {
              warning("qoffset_x is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname qoffset_x-methods
#' @aliases qoffset.x,nifti-method
#' @export
setGeneric("qoffset.x", function(object) standardGeneric("qoffset.x"))
#' @rdname qoffset_x-methods
#' @aliases qoffset.x,nifti-method
#' @export
setMethod("qoffset.x", "nifti", function(object) { object@"qoffset_x" })
#' @rdname qoffset_x-methods
#' @aliases qoffset.x<- 
#' @export
setGeneric("qoffset.x<-", function(object, value) { standardGeneric("qoffset.x<-") })
#' @rdname qoffset_x-methods
#' @aliases qoffset.x<-,nifti-method
#' @export
setMethod("qoffset.x<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "qoffset_x" %in% slotNames(object) ){
              object@"qoffset_x" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("qoffset_x <-", value))               
            } else {
              warning("qoffset_x is not in slotNames of object")
            }                       
            return(object)
          })
