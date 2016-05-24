#' @name intent_p2-methods
#' @title Extract Image Attribute \code{intent_p2}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{intent_p2} field.  
#' @description Methods that act on the \code{intent_p2} field in the
#' NIfTI/ANALYZE header.
#' @rdname intent_p2-methods
#' @aliases intent_p2-methods, intent_p2
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
setGeneric("intent_p2", function(object) standardGeneric("intent_p2"))
#' @rdname intent_p2-methods
#' @aliases intent_p2,nifti-method
#' @export
setMethod("intent_p2", "nifti", function(object) { object@"intent_p2" })
#' @rdname intent_p2-methods
#' @aliases intent_p2<- 
#' @export
setGeneric("intent_p2<-", function(object, value) { standardGeneric("intent_p2<-") })
#' @rdname intent_p2-methods
#' @aliases intent_p2<-,nifti-method
#' @export
setMethod("intent_p2<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "intent_p2" %in% slotNames(object) ){
              object@"intent_p2" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("intent_p2 <-", value))               
            } else {
              warning("intent_p2 is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname intent_p2-methods
#' @aliases intent.p2,nifti-method
#' @export
setGeneric("intent.p2", function(object) standardGeneric("intent.p2"))
#' @rdname intent_p2-methods
#' @aliases intent.p2,nifti-method
#' @export
setMethod("intent.p2", "nifti", function(object) { object@"intent_p2" })
#' @rdname intent_p2-methods
#' @aliases intent.p2<- 
#' @export
setGeneric("intent.p2<-", function(object, value) { standardGeneric("intent.p2<-") })
#' @rdname intent_p2-methods
#' @aliases intent.p2<-,nifti-method
#' @export
setMethod("intent.p2<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "intent_p2" %in% slotNames(object) ){
              object@"intent_p2" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("intent_p2 <-", value))               
            } else {
              warning("intent_p2 is not in slotNames of object")
            }                       
            return(object)
          })
