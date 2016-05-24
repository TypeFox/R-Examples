#' @name intent_code-methods
#' @title Extract Image Attribute \code{intent_code}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{intent_code} field.  
#' @description Methods that act on the \code{intent_code} field in the
#' NIfTI/ANALYZE header.
#' @rdname intent_code-methods
#' @aliases intent_code-methods, intent_code
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
setGeneric("intent_code", function(object) standardGeneric("intent_code"))
#' @rdname intent_code-methods
#' @aliases intent_code,nifti-method
#' @export
setMethod("intent_code", "nifti", function(object) { object@"intent_code" })
#' @rdname intent_code-methods
#' @aliases intent_code<- 
#' @export
setGeneric("intent_code<-", function(object, value) { standardGeneric("intent_code<-") })
#' @rdname intent_code-methods
#' @aliases intent_code<-,nifti-method
#' @export
setMethod("intent_code<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "intent_code" %in% slotNames(object) ){
              object@"intent_code" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("intent_code <-", value))               
            } else {
              warning("intent_code is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname intent_code-methods
#' @aliases intent.code,nifti-method
#' @export
setGeneric("intent.code", function(object) standardGeneric("intent.code"))
#' @rdname intent_code-methods
#' @aliases intent.code,nifti-method
#' @export
setMethod("intent.code", "nifti", function(object) { object@"intent_code" })
#' @rdname intent_code-methods
#' @aliases intent.code<- 
#' @export
setGeneric("intent.code<-", function(object, value) { standardGeneric("intent.code<-") })
#' @rdname intent_code-methods
#' @aliases intent.code<-,nifti-method
#' @export
setMethod("intent.code<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "intent_code" %in% slotNames(object) ){
              object@"intent_code" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("intent_code <-", value))               
            } else {
              warning("intent_code is not in slotNames of object")
            }                       
            return(object)
          })
