#' @name omax-methods
#' @title Extract Image Attribute \code{omax}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{omax} field.  
#' @description Methods that act on the \code{omax} field in the
#' NIfTI/ANALYZE header.
#' @rdname omax-methods
#' @aliases omax-methods, omax
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
setGeneric("omax", function(object) standardGeneric("omax"))
#' @rdname omax-methods
#' @aliases omax,anlz-method
#' @export
setMethod("omax", "anlz", function(object) { object@"omax" })
#' @rdname omax-methods
#' @aliases omax<- 
#' @export
setGeneric("omax<-", function(object, value) { standardGeneric("omax<-") })
#' @rdname omax-methods
#' @aliases omax<-,anlz-method
#' @export
setMethod("omax<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "omax" %in% slotNames(object) ){
              object@"omax" <- value
            } else {
              warning("omax is not in slotNames of object")
            }
            return(object)
          })
