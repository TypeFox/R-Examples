#' @name session_error-methods
#' @title Extract Image Attribute \code{session_error}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{session_error} field.  
#' @description Methods that act on the \code{session_error} field in the
#' NIfTI/ANALYZE header.
#' @rdname session_error-methods
#' @aliases session_error-methods, session_error
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
setGeneric("session_error", function(object) standardGeneric("session_error"))
#' @rdname session_error-methods
#' @aliases session_error,nifti-method
#' @export
setMethod("session_error", "nifti", function(object) { object@"session_error" })
#' @rdname session_error-methods
#' @aliases session_error,anlz-method
#' @export
setMethod("session_error", "anlz", function(object) { object@"session_error" })
#' @rdname session_error-methods
#' @aliases session_error<- 
#' @export
setGeneric("session_error<-", function(object, value) { standardGeneric("session_error<-") })
#' @rdname session_error-methods
#' @aliases session_error<-,nifti-method
#' @export
setMethod("session_error<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "session_error" %in% slotNames(object) ){
              object@"session_error" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("session_error <-", value))               
            } else {
              warning("session_error is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname session_error-methods
#' @aliases session_error<-,anlz-method
#' @export
setMethod("session_error<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "session_error" %in% slotNames(object) ){
              object@"session_error" <- value
            } else {
              warning("session_error is not in slotNames of object")
            }
            return(object)
          })
#' @rdname session_error-methods
#' @aliases session.error,nifti-method
#' @export
setGeneric("session.error", function(object) standardGeneric("session.error"))
#' @rdname session_error-methods
#' @aliases session.error,nifti-method
#' @export
setMethod("session.error", "nifti", function(object) { object@"session_error" })
#' @rdname session_error-methods
#' @aliases session.error,anlz-method
#' @export
setMethod("session.error", "anlz", function(object) { object@"session_error" })
#' @rdname session_error-methods
#' @aliases session.error<- 
#' @export
setGeneric("session.error<-", function(object, value) { standardGeneric("session.error<-") })
#' @rdname session_error-methods
#' @aliases session.error<-,nifti-method
#' @export
setMethod("session.error<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "session_error" %in% slotNames(object) ){
              object@"session_error" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("session_error <-", value))               
            } else {
              warning("session_error is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname session_error-methods
#' @aliases session.error<-,anlz-method
#' @export
setMethod("session.error<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "session_error" %in% slotNames(object) ){
              object@"session_error" <- value
            } else {
              warning("session_error is not in slotNames of object")
            }
            return(object)
          })
