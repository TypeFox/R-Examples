#' @name bitpix-methods
#' @title Extract Image Attribute \code{bitpix}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{bitpix} field.  
#' @description Methods that act on the \code{bitpix} field in the
#' NIfTI/ANALYZE header.
#' @rdname bitpix-methods
#' @aliases bitpix-methods, bitpix
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
setGeneric("bitpix", function(object) standardGeneric("bitpix"))
#' @rdname bitpix-methods
#' @aliases bitpix,nifti-method
#' @export
setMethod("bitpix", "nifti", function(object) { object@"bitpix" })
#' @rdname bitpix-methods
#' @aliases bitpix,anlz-method
#' @export
setMethod("bitpix", "anlz", function(object) { object@"bitpix" })
#' @rdname bitpix-methods
#' @aliases bitpix<- 
#' @export
setGeneric("bitpix<-", function(object, value) { standardGeneric("bitpix<-") })
#' @rdname bitpix-methods
#' @aliases bitpix<-,nifti-method
#' @export
setMethod("bitpix<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "bitpix" %in% slotNames(object) ){
              object@"bitpix" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("bitpix <-", value))               
            } else {
              warning("bitpix is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname bitpix-methods
#' @aliases bitpix<-,anlz-method
#' @export
setMethod("bitpix<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "bitpix" %in% slotNames(object) ){
              object@"bitpix" <- value
            } else {
              warning("bitpix is not in slotNames of object")
            }
            return(object)
          })
