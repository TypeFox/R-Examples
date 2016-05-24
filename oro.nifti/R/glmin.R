#' @name glmin-methods
#' @title Extract Image Attribute \code{glmin}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{glmin} field.  
#' @description Methods that act on the \code{glmin} field in the
#' NIfTI/ANALYZE header.
#' @rdname glmin-methods
#' @aliases glmin-methods, glmin
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
setGeneric("glmin", function(object) standardGeneric("glmin"))
#' @rdname glmin-methods
#' @aliases glmin,nifti-method
#' @export
setMethod("glmin", "nifti", function(object) { object@"glmin" })
#' @rdname glmin-methods
#' @aliases glmin,anlz-method
#' @export
setMethod("glmin", "anlz", function(object) { object@"glmin" })
#' @rdname glmin-methods
#' @aliases glmin<- 
#' @export
setGeneric("glmin<-", function(object, value) { standardGeneric("glmin<-") })
#' @rdname glmin-methods
#' @aliases glmin<-,nifti-method
#' @export
setMethod("glmin<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "glmin" %in% slotNames(object) ){
              object@"glmin" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("glmin <-", value))               
            } else {
              warning("glmin is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname glmin-methods
#' @aliases glmin<-,anlz-method
#' @export
setMethod("glmin<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "glmin" %in% slotNames(object) ){
              object@"glmin" <- value
            } else {
              warning("glmin is not in slotNames of object")
            }
            return(object)
          })
