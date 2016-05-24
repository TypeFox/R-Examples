#' @name quatern_c-methods
#' @title Extract Image Attribute \code{quatern_c}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{quatern_c} field.  
#' @description Methods that act on the \code{quatern_c} field in the
#' NIfTI/ANALYZE header.
#' @rdname quatern_c-methods
#' @aliases quatern_c-methods, quatern_c
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
setGeneric("quatern_c", function(object) standardGeneric("quatern_c"))
#' @rdname quatern_c-methods
#' @aliases quatern_c,nifti-method
#' @export
setMethod("quatern_c", "nifti", function(object) { object@"quatern_c" })
#' @rdname quatern_c-methods
#' @aliases quatern_c<- 
#' @export
setGeneric("quatern_c<-", function(object, value) { standardGeneric("quatern_c<-") })
#' @rdname quatern_c-methods
#' @aliases quatern_c<-,nifti-method
#' @export
setMethod("quatern_c<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "quatern_c" %in% slotNames(object) ){
              object@"quatern_c" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("quatern_c <-", value))               
            } else {
              warning("quatern_c is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname quatern_c-methods
#' @aliases quatern.c,nifti-method
#' @export
setGeneric("quatern.c", function(object) standardGeneric("quatern.c"))
#' @rdname quatern_c-methods
#' @aliases quatern.c,nifti-method
#' @export
setMethod("quatern.c", "nifti", function(object) { object@"quatern_c" })
#' @rdname quatern_c-methods
#' @aliases quatern.c<- 
#' @export
setGeneric("quatern.c<-", function(object, value) { standardGeneric("quatern.c<-") })
#' @rdname quatern_c-methods
#' @aliases quatern.c<-,nifti-method
#' @export
setMethod("quatern.c<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "quatern_c" %in% slotNames(object) ){
              object@"quatern_c" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("quatern_c <-", value))               
            } else {
              warning("quatern_c is not in slotNames of object")
            }                       
            return(object)
          })
