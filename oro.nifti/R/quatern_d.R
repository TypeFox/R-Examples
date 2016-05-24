#' @name quatern_d-methods
#' @title Extract Image Attribute \code{quatern_d}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{quatern_d} field.  
#' @description Methods that act on the \code{quatern_d} field in the
#' NIfTI/ANALYZE header.
#' @rdname quatern_d-methods
#' @aliases quatern_d-methods, quatern_d
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
setGeneric("quatern_d", function(object) standardGeneric("quatern_d"))
#' @rdname quatern_d-methods
#' @aliases quatern_d,nifti-method
#' @export
setMethod("quatern_d", "nifti", function(object) { object@"quatern_d" })
#' @rdname quatern_d-methods
#' @aliases quatern_d<- 
#' @export
setGeneric("quatern_d<-", function(object, value) { standardGeneric("quatern_d<-") })
#' @rdname quatern_d-methods
#' @aliases quatern_d<-,nifti-method
#' @export
setMethod("quatern_d<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "quatern_d" %in% slotNames(object) ){
              object@"quatern_d" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("quatern_d <-", value))               
            } else {
              warning("quatern_d is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname quatern_d-methods
#' @aliases quatern.d,nifti-method
#' @export
setGeneric("quatern.d", function(object) standardGeneric("quatern.d"))
#' @rdname quatern_d-methods
#' @aliases quatern.d,nifti-method
#' @export
setMethod("quatern.d", "nifti", function(object) { object@"quatern_d" })
#' @rdname quatern_d-methods
#' @aliases quatern.d<- 
#' @export
setGeneric("quatern.d<-", function(object, value) { standardGeneric("quatern.d<-") })
#' @rdname quatern_d-methods
#' @aliases quatern.d<-,nifti-method
#' @export
setMethod("quatern.d<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "quatern_d" %in% slotNames(object) ){
              object@"quatern_d" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("quatern_d <-", value))               
            } else {
              warning("quatern_d is not in slotNames of object")
            }                       
            return(object)
          })
