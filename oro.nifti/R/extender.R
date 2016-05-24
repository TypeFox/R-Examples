#' @name extender-methods
#' @title Extract Image Attribute \code{extender}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{extender} field.  
#' @description Methods that act on the \code{extender} field in the
#' NIfTI/ANALYZE header.
#' @rdname extender-methods
#' @aliases extender-methods, extender
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
setGeneric("extender", function(object) standardGeneric("extender"))
#' @rdname extender-methods
#' @aliases extender,nifti-method
#' @export
setMethod("extender", "nifti", function(object) { object@"extender" })
#' @rdname extender-methods
#' @aliases extender<- 
#' @export
setGeneric("extender<-", function(object, value) { standardGeneric("extender<-") })
#' @rdname extender-methods
#' @aliases extender<-,nifti-method
#' @export
setMethod("extender<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "extender" %in% slotNames(object) ){
              object@"extender" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("extender <-", value))               
            } else {
              warning("extender is not in slotNames of object")
            }                       
            return(object)
          })
