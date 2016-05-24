#' @name slice_code-methods
#' @title Extract Image Attribute \code{slice_code}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{slice_code} field.  
#' @description Methods that act on the \code{slice_code} field in the
#' NIfTI/ANALYZE header.
#' @rdname slice_code-methods
#' @aliases slice_code-methods, slice_code
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
setGeneric("slice_code", function(object) standardGeneric("slice_code"))
#' @rdname slice_code-methods
#' @aliases slice_code,nifti-method
#' @export
setMethod("slice_code", "nifti", function(object) { object@"slice_code" })
#' @rdname slice_code-methods
#' @aliases slice_code<- 
#' @export
setGeneric("slice_code<-", function(object, value) { standardGeneric("slice_code<-") })
#' @rdname slice_code-methods
#' @aliases slice_code<-,nifti-method
#' @export
setMethod("slice_code<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "slice_code" %in% slotNames(object) ){
              object@"slice_code" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("slice_code <-", value))               
            } else {
              warning("slice_code is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname slice_code-methods
#' @aliases slice.code,nifti-method
#' @export
setGeneric("slice.code", function(object) standardGeneric("slice.code"))
#' @rdname slice_code-methods
#' @aliases slice.code,nifti-method
#' @export
setMethod("slice.code", "nifti", function(object) { object@"slice_code" })
#' @rdname slice_code-methods
#' @aliases slice.code<- 
#' @export
setGeneric("slice.code<-", function(object, value) { standardGeneric("slice.code<-") })
#' @rdname slice_code-methods
#' @aliases slice.code<-,nifti-method
#' @export
setMethod("slice.code<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "slice_code" %in% slotNames(object) ){
              object@"slice_code" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("slice_code <-", value))               
            } else {
              warning("slice_code is not in slotNames of object")
            }                       
            return(object)
          })
