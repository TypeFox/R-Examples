#' @name funused1-methods
#' @title Extract Image Attribute \code{funused1}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{funused1} field.  
#' @description Methods that act on the \code{funused1} field in the
#' NIfTI/ANALYZE header.
#' @rdname funused1-methods
#' @aliases funused1-methods, funused1
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
setGeneric("funused1", function(object) standardGeneric("funused1"))
#' @rdname funused1-methods
#' @aliases funused1,anlz-method
#' @export
setMethod("funused1", "anlz", function(object) { object@"funused1" })
#' @rdname funused1-methods
#' @aliases funused1<- 
#' @export
setGeneric("funused1<-", function(object, value) { standardGeneric("funused1<-") })
#' @rdname funused1-methods
#' @aliases funused1<-,anlz-method
#' @export
setMethod("funused1<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "funused1" %in% slotNames(object) ){
              object@"funused1" <- value
            } else {
              warning("funused1 is not in slotNames of object")
            }
            return(object)
          })
