#' @name intent_name-methods
#' @title Extract Image Attribute \code{intent_name}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{intent_name} field.  
#' @description Methods that act on the \code{intent_name} field in the
#' NIfTI/ANALYZE header.
#' @rdname intent_name-methods
#' @aliases intent_name-methods, intent_name
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
setGeneric("intent_name", function(object) standardGeneric("intent_name"))
#' @rdname intent_name-methods
#' @aliases intent_name,nifti-method
#' @export
setMethod("intent_name", "nifti", function(object) { object@"intent_name" })
#' @rdname intent_name-methods
#' @aliases intent_name<- 
#' @export
setGeneric("intent_name<-", function(object, value) { standardGeneric("intent_name<-") })
#' @rdname intent_name-methods
#' @aliases intent_name<-,nifti-method
#' @export
setMethod("intent_name<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "intent_name" %in% slotNames(object) ){
              object@"intent_name" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("intent_name <-", value))               
            } else {
              warning("intent_name is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname intent_name-methods
#' @aliases intent.name,nifti-method
#' @export
setGeneric("intent.name", function(object) standardGeneric("intent.name"))
#' @rdname intent_name-methods
#' @aliases intent.name,nifti-method
#' @export
setMethod("intent.name", "nifti", function(object) { object@"intent_name" })
#' @rdname intent_name-methods
#' @aliases intent.name<- 
#' @export
setGeneric("intent.name<-", function(object, value) { standardGeneric("intent.name<-") })
#' @rdname intent_name-methods
#' @aliases intent.name<-,nifti-method
#' @export
setMethod("intent.name<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "intent_name" %in% slotNames(object) ){
              object@"intent_name" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("intent_name <-", value))               
            } else {
              warning("intent_name is not in slotNames of object")
            }                       
            return(object)
          })
