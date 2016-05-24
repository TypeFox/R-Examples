#' @name slice_start-methods
#' @title Extract Image Attribute \code{slice_start}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{slice_start} field.  
#' @description Methods that act on the \code{slice_start} field in the
#' NIfTI/ANALYZE header.
#' @rdname slice_start-methods
#' @aliases slice_start-methods, slice_start
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
setGeneric("slice_start", function(object) standardGeneric("slice_start"))
#' @rdname slice_start-methods
#' @aliases slice_start,nifti-method
#' @export
setMethod("slice_start", "nifti", function(object) { object@"slice_start" })
#' @rdname slice_start-methods
#' @aliases slice_start<- 
#' @export
setGeneric("slice_start<-", function(object, value) { standardGeneric("slice_start<-") })
#' @rdname slice_start-methods
#' @aliases slice_start<-,nifti-method
#' @export
setMethod("slice_start<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "slice_start" %in% slotNames(object) ){
              object@"slice_start" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("slice_start <-", value))               
            } else {
              warning("slice_start is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname slice_start-methods
#' @aliases slice.start,nifti-method
#' @export
setGeneric("slice.start", function(object) standardGeneric("slice.start"))
#' @rdname slice_start-methods
#' @aliases slice.start,nifti-method
#' @export
setMethod("slice.start", "nifti", function(object) { object@"slice_start" })
#' @rdname slice_start-methods
#' @aliases slice.start<- 
#' @export
setGeneric("slice.start<-", function(object, value) { standardGeneric("slice.start<-") })
#' @rdname slice_start-methods
#' @aliases slice.start<-,nifti-method
#' @export
setMethod("slice.start<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "slice_start" %in% slotNames(object) ){
              object@"slice_start" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("slice_start <-", value))               
            } else {
              warning("slice_start is not in slotNames of object")
            }                       
            return(object)
          })
