#' @name srow_x-methods
#' @title Extract Image Attribute \code{srow_x}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{srow_x} field.  
#' @description Methods that act on the \code{srow_x} field in the
#' NIfTI/ANALYZE header.
#' @rdname srow_x-methods
#' @aliases srow_x-methods, srow_x
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
setGeneric("srow_x", function(object) standardGeneric("srow_x"))
#' @rdname srow_x-methods
#' @aliases srow_x,nifti-method
#' @export
setMethod("srow_x", "nifti", function(object) { object@"srow_x" })
#' @rdname srow_x-methods
#' @aliases srow_x<- 
#' @export
setGeneric("srow_x<-", function(object, value) { standardGeneric("srow_x<-") })
#' @rdname srow_x-methods
#' @aliases srow_x<-,nifti-method
#' @export
setMethod("srow_x<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "srow_x" %in% slotNames(object) ){
              object@"srow_x" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("srow_x <-", value))               
            } else {
              warning("srow_x is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname srow_x-methods
#' @aliases srow.x,nifti-method
#' @export
setGeneric("srow.x", function(object) standardGeneric("srow.x"))
#' @rdname srow_x-methods
#' @aliases srow.x,nifti-method
#' @export
setMethod("srow.x", "nifti", function(object) { object@"srow_x" })
#' @rdname srow_x-methods
#' @aliases srow.x<- 
#' @export
setGeneric("srow.x<-", function(object, value) { standardGeneric("srow.x<-") })
#' @rdname srow_x-methods
#' @aliases srow.x<-,nifti-method
#' @export
setMethod("srow.x<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "srow_x" %in% slotNames(object) ){
              object@"srow_x" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("srow_x <-", value))               
            } else {
              warning("srow_x is not in slotNames of object")
            }                       
            return(object)
          })
