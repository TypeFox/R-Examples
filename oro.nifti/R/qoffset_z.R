#' @name qoffset_z-methods
#' @title Extract Image Attribute \code{qoffset_z}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{qoffset_z} field.  
#' @description Methods that act on the \code{qoffset_z} field in the
#' NIfTI/ANALYZE header.
#' @rdname qoffset_z-methods
#' @aliases qoffset_z-methods, qoffset_z
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
setGeneric("qoffset_z", function(object) standardGeneric("qoffset_z"))
#' @rdname qoffset_z-methods
#' @aliases qoffset_z,nifti-method
#' @export
setMethod("qoffset_z", "nifti", function(object) { object@"qoffset_z" })
#' @rdname qoffset_z-methods
#' @aliases qoffset_z<- 
#' @export
setGeneric("qoffset_z<-", function(object, value) { standardGeneric("qoffset_z<-") })
#' @rdname qoffset_z-methods
#' @aliases qoffset_z<-,nifti-method
#' @export
setMethod("qoffset_z<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "qoffset_z" %in% slotNames(object) ){
              object@"qoffset_z" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("qoffset_z <-", value))               
            } else {
              warning("qoffset_z is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname qoffset_z-methods
#' @aliases qoffset.z,nifti-method
#' @export
setGeneric("qoffset.z", function(object) standardGeneric("qoffset.z"))
#' @rdname qoffset_z-methods
#' @aliases qoffset.z,nifti-method
#' @export
setMethod("qoffset.z", "nifti", function(object) { object@"qoffset_z" })
#' @rdname qoffset_z-methods
#' @aliases qoffset.z<- 
#' @export
setGeneric("qoffset.z<-", function(object, value) { standardGeneric("qoffset.z<-") })
#' @rdname qoffset_z-methods
#' @aliases qoffset.z<-,nifti-method
#' @export
setMethod("qoffset.z<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "qoffset_z" %in% slotNames(object) ){
              object@"qoffset_z" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("qoffset_z <-", value))               
            } else {
              warning("qoffset_z is not in slotNames of object")
            }                       
            return(object)
          })
