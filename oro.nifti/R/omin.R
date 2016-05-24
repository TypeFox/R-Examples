#' @name omin-methods
#' @title Extract Image Attribute \code{omin}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{omin} field.  
#' @description Methods that act on the \code{omin} field in the
#' NIfTI/ANALYZE header.
#' @rdname omin-methods
#' @aliases omin-methods, omin
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
setGeneric("omin", function(object) standardGeneric("omin"))
#' @rdname omin-methods
#' @aliases omin,anlz-method
#' @export
setMethod("omin", "anlz", function(object) { object@"omin" })
#' @rdname omin-methods
#' @aliases omin<- 
#' @export
setGeneric("omin<-", function(object, value) { standardGeneric("omin<-") })
#' @rdname omin-methods
#' @aliases omin<-,anlz-method
#' @export
setMethod("omin<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "omin" %in% slotNames(object) ){
              object@"omin" <- value
            } else {
              warning("omin is not in slotNames of object")
            }
            return(object)
          })
