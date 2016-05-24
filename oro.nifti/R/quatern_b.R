#' @name quatern_b-methods
#' @title Extract Image Attribute \code{quatern_b}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{quatern_b} field.  
#' @description Methods that act on the \code{quatern_b} field in the
#' NIfTI/ANALYZE header.
#' @rdname quatern_b-methods
#' @aliases quatern_b-methods, quatern_b
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
setGeneric("quatern_b", function(object) standardGeneric("quatern_b"))
#' @rdname quatern_b-methods
#' @aliases quatern_b,nifti-method
#' @export
setMethod("quatern_b", "nifti", function(object) { object@"quatern_b" })
#' @rdname quatern_b-methods
#' @aliases quatern_b<- 
#' @export
setGeneric("quatern_b<-", function(object, value) { standardGeneric("quatern_b<-") })
#' @rdname quatern_b-methods
#' @aliases quatern_b<-,nifti-method
#' @export
setMethod("quatern_b<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "quatern_b" %in% slotNames(object) ){
              object@"quatern_b" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("quatern_b <-", value))               
            } else {
              warning("quatern_b is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname quatern_b-methods
#' @aliases quatern.b,nifti-method
#' @export
setGeneric("quatern.b", function(object) standardGeneric("quatern.b"))
#' @rdname quatern_b-methods
#' @aliases quatern.b,nifti-method
#' @export
setMethod("quatern.b", "nifti", function(object) { object@"quatern_b" })
#' @rdname quatern_b-methods
#' @aliases quatern.b<- 
#' @export
setGeneric("quatern.b<-", function(object, value) { standardGeneric("quatern.b<-") })
#' @rdname quatern_b-methods
#' @aliases quatern.b<-,nifti-method
#' @export
setMethod("quatern.b<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "quatern_b" %in% slotNames(object) ){
              object@"quatern_b" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("quatern_b <-", value))               
            } else {
              warning("quatern_b is not in slotNames of object")
            }                       
            return(object)
          })
