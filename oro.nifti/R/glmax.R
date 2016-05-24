#' @name glmax-methods
#' @title Extract Image Attribute \code{glmax}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{glmax} field.  
#' @description Methods that act on the \code{glmax} field in the
#' NIfTI/ANALYZE header.
#' @rdname glmax-methods
#' @aliases glmax-methods, glmax
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
setGeneric("glmax", function(object) standardGeneric("glmax"))
#' @rdname glmax-methods
#' @aliases glmax,nifti-method
#' @export
setMethod("glmax", "nifti", function(object) { object@"glmax" })
#' @rdname glmax-methods
#' @aliases glmax,anlz-method
#' @export
setMethod("glmax", "anlz", function(object) { object@"glmax" })
#' @rdname glmax-methods
#' @aliases glmax<- 
#' @export
setGeneric("glmax<-", function(object, value) { standardGeneric("glmax<-") })
#' @rdname glmax-methods
#' @aliases glmax<-,nifti-method
#' @export
setMethod("glmax<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "glmax" %in% slotNames(object) ){
              object@"glmax" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("glmax <-", value))               
            } else {
              warning("glmax is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname glmax-methods
#' @aliases glmax<-,anlz-method
#' @export
setMethod("glmax<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "glmax" %in% slotNames(object) ){
              object@"glmax" <- value
            } else {
              warning("glmax is not in slotNames of object")
            }
            return(object)
          })
