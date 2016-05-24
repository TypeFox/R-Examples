#' @name slice_duration-methods
#' @title Extract Image Attribute \code{slice_duration}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{slice_duration} field.  
#' @description Methods that act on the \code{slice_duration} field in the
#' NIfTI/ANALYZE header.
#' @rdname slice_duration-methods
#' @aliases slice_duration-methods, slice_duration
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
setGeneric("slice_duration", function(object) standardGeneric("slice_duration"))
#' @rdname slice_duration-methods
#' @aliases slice_duration,nifti-method
#' @export
setMethod("slice_duration", "nifti", function(object) { object@"slice_duration" })
#' @rdname slice_duration-methods
#' @aliases slice_duration<- 
#' @export
setGeneric("slice_duration<-", function(object, value) { standardGeneric("slice_duration<-") })
#' @rdname slice_duration-methods
#' @aliases slice_duration<-,nifti-method
#' @export
setMethod("slice_duration<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "slice_duration" %in% slotNames(object) ){
              object@"slice_duration" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("slice_duration <-", value))               
            } else {
              warning("slice_duration is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname slice_duration-methods
#' @aliases slice.duration,nifti-method
#' @export
setGeneric("slice.duration", function(object) standardGeneric("slice.duration"))
#' @rdname slice_duration-methods
#' @aliases slice.duration,nifti-method
#' @export
setMethod("slice.duration", "nifti", function(object) { object@"slice_duration" })
#' @rdname slice_duration-methods
#' @aliases slice.duration<- 
#' @export
setGeneric("slice.duration<-", function(object, value) { standardGeneric("slice.duration<-") })
#' @rdname slice_duration-methods
#' @aliases slice.duration<-,nifti-method
#' @export
setMethod("slice.duration<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "slice_duration" %in% slotNames(object) ){
              object@"slice_duration" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("slice_duration <-", value))               
            } else {
              warning("slice_duration is not in slotNames of object")
            }                       
            return(object)
          })
