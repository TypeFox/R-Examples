#' @name srow_z-methods
#' @title Extract Image Attribute \code{srow_z}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{srow_z} field.  
#' @description Methods that act on the \code{srow_z} field in the
#' NIfTI/ANALYZE header.
#' @rdname srow_z-methods
#' @aliases srow_z-methods, srow_z
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
setGeneric("srow_z", function(object) standardGeneric("srow_z"))
#' @rdname srow_z-methods
#' @aliases srow_z,nifti-method
#' @export
setMethod("srow_z", "nifti", function(object) { object@"srow_z" })
#' @rdname srow_z-methods
#' @aliases srow_z<- 
#' @export
setGeneric("srow_z<-", function(object, value) { standardGeneric("srow_z<-") })
#' @rdname srow_z-methods
#' @aliases srow_z<-,nifti-method
#' @export
setMethod("srow_z<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "srow_z" %in% slotNames(object) ){
              object@"srow_z" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("srow_z <-", value))               
            } else {
              warning("srow_z is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname srow_z-methods
#' @aliases srow.z,nifti-method
#' @export
setGeneric("srow.z", function(object) standardGeneric("srow.z"))
#' @rdname srow_z-methods
#' @aliases srow.z,nifti-method
#' @export
setMethod("srow.z", "nifti", function(object) { object@"srow_z" })
#' @rdname srow_z-methods
#' @aliases srow.z<- 
#' @export
setGeneric("srow.z<-", function(object, value) { standardGeneric("srow.z<-") })
#' @rdname srow_z-methods
#' @aliases srow.z<-,nifti-method
#' @export
setMethod("srow.z<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "srow_z" %in% slotNames(object) ){
              object@"srow_z" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("srow_z <-", value))               
            } else {
              warning("srow_z is not in slotNames of object")
            }                       
            return(object)
          })
