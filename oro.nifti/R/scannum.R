#' @name scannum-methods
#' @title Extract Image Attribute \code{scannum}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{scannum} field.  
#' @description Methods that act on the \code{scannum} field in the
#' NIfTI/ANALYZE header.
#' @rdname scannum-methods
#' @aliases scannum-methods, scannum
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
setGeneric("scannum", function(object) standardGeneric("scannum"))
#' @rdname scannum-methods
#' @aliases scannum,anlz-method
#' @export
setMethod("scannum", "anlz", function(object) { object@"scannum" })
#' @rdname scannum-methods
#' @aliases scannum<- 
#' @export
setGeneric("scannum<-", function(object, value) { standardGeneric("scannum<-") })
#' @rdname scannum-methods
#' @aliases scannum<-,anlz-method
#' @export
setMethod("scannum<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "scannum" %in% slotNames(object) ){
              object@"scannum" <- value
            } else {
              warning("scannum is not in slotNames of object")
            }
            return(object)
          })
