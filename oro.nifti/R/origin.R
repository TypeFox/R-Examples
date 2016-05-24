#' @name origin-methods
#' @title Extract Image Attribute \code{origin}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{origin} field.  
#' @description Methods that act on the \code{origin} field in the
#' NIfTI/ANALYZE header.
#' @rdname origin-methods
#' @aliases origin-methods, origin
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
setGeneric("origin", function(object) standardGeneric("origin"))
#' @rdname origin-methods
#' @aliases origin,anlz-method
#' @export
setMethod("origin", "anlz", function(object) { object@"origin" })
#' @rdname origin-methods
#' @aliases origin<- 
#' @export
setGeneric("origin<-", function(object, value) { standardGeneric("origin<-") })
#' @rdname origin-methods
#' @aliases origin<-,anlz-method
#' @export
setMethod("origin<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "origin" %in% slotNames(object) ){
              object@"origin" <- value
            } else {
              warning("origin is not in slotNames of object")
            }
            return(object)
          })
