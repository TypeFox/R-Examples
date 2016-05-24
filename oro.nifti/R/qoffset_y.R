#' @name qoffset_y-methods
#' @title Extract Image Attribute \code{qoffset_y}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{qoffset_y} field.  
#' @description Methods that act on the \code{qoffset_y} field in the
#' NIfTI/ANALYZE header.
#' @rdname qoffset_y-methods
#' @aliases qoffset_y-methods, qoffset_y
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
setGeneric("qoffset_y", function(object) standardGeneric("qoffset_y"))
#' @rdname qoffset_y-methods
#' @aliases qoffset_y,nifti-method
#' @export
setMethod("qoffset_y", "nifti", function(object) { object@"qoffset_y" })
#' @rdname qoffset_y-methods
#' @aliases qoffset_y<- 
#' @export
setGeneric("qoffset_y<-", function(object, value) { standardGeneric("qoffset_y<-") })
#' @rdname qoffset_y-methods
#' @aliases qoffset_y<-,nifti-method
#' @export
setMethod("qoffset_y<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "qoffset_y" %in% slotNames(object) ){
              object@"qoffset_y" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("qoffset_y <-", value))               
            } else {
              warning("qoffset_y is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname qoffset_y-methods
#' @aliases qoffset.y,nifti-method
#' @export
setGeneric("qoffset.y", function(object) standardGeneric("qoffset.y"))
#' @rdname qoffset_y-methods
#' @aliases qoffset.y,nifti-method
#' @export
setMethod("qoffset.y", "nifti", function(object) { object@"qoffset_y" })
#' @rdname qoffset_y-methods
#' @aliases qoffset.y<- 
#' @export
setGeneric("qoffset.y<-", function(object, value) { standardGeneric("qoffset.y<-") })
#' @rdname qoffset_y-methods
#' @aliases qoffset.y<-,nifti-method
#' @export
setMethod("qoffset.y<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "qoffset_y" %in% slotNames(object) ){
              object@"qoffset_y" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("qoffset_y <-", value))               
            } else {
              warning("qoffset_y is not in slotNames of object")
            }                       
            return(object)
          })
