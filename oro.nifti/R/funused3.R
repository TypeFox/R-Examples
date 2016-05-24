#' @name funused3-methods
#' @title Extract Image Attribute \code{funused3}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{funused3} field.  
#' @description Methods that act on the \code{funused3} field in the
#' NIfTI/ANALYZE header.
#' @rdname funused3-methods
#' @aliases funused3-methods, funused3
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
setGeneric("funused3", function(object) standardGeneric("funused3"))
#' @rdname funused3-methods
#' @aliases funused3,anlz-method
#' @export
setMethod("funused3", "anlz", function(object) { object@"funused3" })
#' @rdname funused3-methods
#' @aliases funused3<- 
#' @export
setGeneric("funused3<-", function(object, value) { standardGeneric("funused3<-") })
#' @rdname funused3-methods
#' @aliases funused3<-,anlz-method
#' @export
setMethod("funused3<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "funused3" %in% slotNames(object) ){
              object@"funused3" <- value
            } else {
              warning("funused3 is not in slotNames of object")
            }
            return(object)
          })
