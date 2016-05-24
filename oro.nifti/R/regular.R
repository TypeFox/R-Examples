#' @name regular-methods
#' @title Extract Image Attribute \code{regular}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{regular} field.  
#' @description Methods that act on the \code{regular} field in the
#' NIfTI/ANALYZE header.
#' @rdname regular-methods
#' @aliases regular-methods, regular
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
setGeneric("regular", function(object) standardGeneric("regular"))
#' @rdname regular-methods
#' @aliases regular,nifti-method
#' @export
setMethod("regular", "nifti", function(object) { object@"regular" })
#' @rdname regular-methods
#' @aliases regular,anlz-method
#' @export
setMethod("regular", "anlz", function(object) { object@"regular" })
#' @rdname regular-methods
#' @aliases regular<- 
#' @export
setGeneric("regular<-", function(object, value) { standardGeneric("regular<-") })
#' @rdname regular-methods
#' @aliases regular<-,nifti-method
#' @export
setMethod("regular<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "regular" %in% slotNames(object) ){
              object@"regular" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("regular <-", value))               
            } else {
              warning("regular is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname regular-methods
#' @aliases regular<-,anlz-method
#' @export
setMethod("regular<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "regular" %in% slotNames(object) ){
              object@"regular" <- value
            } else {
              warning("regular is not in slotNames of object")
            }
            return(object)
          })
