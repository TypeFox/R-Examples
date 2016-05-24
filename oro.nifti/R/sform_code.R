#' @name sform_code-methods
#' @title Extract Image Attribute \code{sform_code}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{sform_code} field.  
#' @description Methods that act on the \code{sform_code} field in the
#' NIfTI/ANALYZE header.
#' @rdname sform_code-methods
#' @aliases sform_code-methods, sform_code
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
setGeneric("sform_code", function(object) standardGeneric("sform_code"))
#' @rdname sform_code-methods
#' @aliases sform_code,nifti-method
#' @export
setMethod("sform_code", "nifti", function(object) { object@"sform_code" })
#' @rdname sform_code-methods
#' @aliases sform_code<- 
#' @export
setGeneric("sform_code<-", function(object, value) { standardGeneric("sform_code<-") })
#' @rdname sform_code-methods
#' @aliases sform_code<-,nifti-method
#' @export
setMethod("sform_code<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "sform_code" %in% slotNames(object) ){
              object@"sform_code" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("sform_code <-", value))               
            } else {
              warning("sform_code is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname sform_code-methods
#' @aliases sform.code,nifti-method
#' @export
setGeneric("sform.code", function(object) standardGeneric("sform.code"))
#' @rdname sform_code-methods
#' @aliases sform.code,nifti-method
#' @export
setMethod("sform.code", "nifti", function(object) { object@"sform_code" })
#' @rdname sform_code-methods
#' @aliases sform.code<- 
#' @export
setGeneric("sform.code<-", function(object, value) { standardGeneric("sform.code<-") })
#' @rdname sform_code-methods
#' @aliases sform.code<-,nifti-method
#' @export
setMethod("sform.code<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "sform_code" %in% slotNames(object) ){
              object@"sform_code" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("sform_code <-", value))               
            } else {
              warning("sform_code is not in slotNames of object")
            }                       
            return(object)
          })
