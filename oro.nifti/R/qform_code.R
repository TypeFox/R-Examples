#' @name qform_code-methods
#' @title Extract Image Attribute \code{qform_code}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{qform_code} field.  
#' @description Methods that act on the \code{qform_code} field in the
#' NIfTI/ANALYZE header.
#' @rdname qform_code-methods
#' @aliases qform_code-methods, qform_code
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
setGeneric("qform_code", function(object) standardGeneric("qform_code"))
#' @rdname qform_code-methods
#' @aliases qform_code,nifti-method
#' @export
setMethod("qform_code", "nifti", function(object) { object@"qform_code" })
#' @rdname qform_code-methods
#' @aliases qform_code<- 
#' @export
setGeneric("qform_code<-", function(object, value) { standardGeneric("qform_code<-") })
#' @rdname qform_code-methods
#' @aliases qform_code<-,nifti-method
#' @export
setMethod("qform_code<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "qform_code" %in% slotNames(object) ){
              object@"qform_code" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("qform_code <-", value))               
            } else {
              warning("qform_code is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname qform_code-methods
#' @aliases qform.code,nifti-method
#' @export
setGeneric("qform.code", function(object) standardGeneric("qform.code"))
#' @rdname qform_code-methods
#' @aliases qform.code,nifti-method
#' @export
setMethod("qform.code", "nifti", function(object) { object@"qform_code" })
#' @rdname qform_code-methods
#' @aliases qform.code<- 
#' @export
setGeneric("qform.code<-", function(object, value) { standardGeneric("qform.code<-") })
#' @rdname qform_code-methods
#' @aliases qform.code<-,nifti-method
#' @export
setMethod("qform.code<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "qform_code" %in% slotNames(object) ){
              object@"qform_code" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("qform_code <-", value))               
            } else {
              warning("qform_code is not in slotNames of object")
            }                       
            return(object)
          })
