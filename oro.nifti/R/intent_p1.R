#' @name intent_p1-methods
#' @title Extract Image Attribute \code{intent_p1}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{intent_p1} field.  
#' @description Methods that act on the \code{intent_p1} field in the
#' NIfTI/ANALYZE header.
#' @rdname intent_p1-methods
#' @aliases intent_p1-methods, intent_p1
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
setGeneric("intent_p1", function(object) standardGeneric("intent_p1"))
#' @rdname intent_p1-methods
#' @aliases intent_p1,nifti-method
#' @export
setMethod("intent_p1", "nifti", function(object) { object@"intent_p1" })
#' @rdname intent_p1-methods
#' @aliases intent_p1<- 
#' @export
setGeneric("intent_p1<-", function(object, value) { standardGeneric("intent_p1<-") })
#' @rdname intent_p1-methods
#' @aliases intent_p1<-,nifti-method
#' @export
setMethod("intent_p1<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "intent_p1" %in% slotNames(object) ){
              object@"intent_p1" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("intent_p1 <-", value))               
            } else {
              warning("intent_p1 is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname intent_p1-methods
#' @aliases intent.p1,nifti-method
#' @export
setGeneric("intent.p1", function(object) standardGeneric("intent.p1"))
#' @rdname intent_p1-methods
#' @aliases intent.p1,nifti-method
#' @export
setMethod("intent.p1", "nifti", function(object) { object@"intent_p1" })
#' @rdname intent_p1-methods
#' @aliases intent.p1<- 
#' @export
setGeneric("intent.p1<-", function(object, value) { standardGeneric("intent.p1<-") })
#' @rdname intent_p1-methods
#' @aliases intent.p1<-,nifti-method
#' @export
setMethod("intent.p1<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "intent_p1" %in% slotNames(object) ){
              object@"intent_p1" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("intent_p1 <-", value))               
            } else {
              warning("intent_p1 is not in slotNames of object")
            }                       
            return(object)
          })
