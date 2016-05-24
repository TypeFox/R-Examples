#' @name extents-methods
#' @title Extract Image Attribute \code{extents}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{extents} field.  
#' @description Methods that act on the \code{extents} field in the
#' NIfTI/ANALYZE header.
#' @rdname extents-methods
#' @aliases extents-methods, extents
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
setGeneric("extents", function(object) standardGeneric("extents"))
#' @rdname extents-methods
#' @aliases extents,nifti-method
#' @export
setMethod("extents", "nifti", function(object) { object@"extents" })
#' @rdname extents-methods
#' @aliases extents,anlz-method
#' @export
setMethod("extents", "anlz", function(object) { object@"extents" })
#' @rdname extents-methods
#' @aliases extents<- 
#' @export
setGeneric("extents<-", function(object, value) { standardGeneric("extents<-") })
#' @rdname extents-methods
#' @aliases extents<-,nifti-method
#' @export
setMethod("extents<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "extents" %in% slotNames(object) ){
              object@"extents" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("extents <-", value))               
            } else {
              warning("extents is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname extents-methods
#' @aliases extents<-,anlz-method
#' @export
setMethod("extents<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "extents" %in% slotNames(object) ){
              object@"extents" <- value
            } else {
              warning("extents is not in slotNames of object")
            }
            return(object)
          })
