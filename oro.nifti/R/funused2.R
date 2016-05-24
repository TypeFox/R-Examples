#' @name funused2-methods
#' @title Extract Image Attribute \code{funused2}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{funused2} field.  
#' @description Methods that act on the \code{funused2} field in the
#' NIfTI/ANALYZE header.
#' @rdname funused2-methods
#' @aliases funused2-methods, funused2
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
setGeneric("funused2", function(object) standardGeneric("funused2"))
#' @rdname funused2-methods
#' @aliases funused2,anlz-method
#' @export
setMethod("funused2", "anlz", function(object) { object@"funused2" })
#' @rdname funused2-methods
#' @aliases funused2<- 
#' @export
setGeneric("funused2<-", function(object, value) { standardGeneric("funused2<-") })
#' @rdname funused2-methods
#' @aliases funused2<-,anlz-method
#' @export
setMethod("funused2<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "funused2" %in% slotNames(object) ){
              object@"funused2" <- value
            } else {
              warning("funused2 is not in slotNames of object")
            }
            return(object)
          })
