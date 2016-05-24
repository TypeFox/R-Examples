#' @name verified-methods
#' @title Extract Image Attribute \code{verified}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{verified} field.  
#' @description Methods that act on the \code{verified} field in the
#' NIfTI/ANALYZE header.
#' @rdname verified-methods
#' @aliases verified-methods, verified
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
setGeneric("verified", function(object) standardGeneric("verified"))
#' @rdname verified-methods
#' @aliases verified,anlz-method
#' @export
setMethod("verified", "anlz", function(object) { object@"verified" })
#' @rdname verified-methods
#' @aliases verified<- 
#' @export
setGeneric("verified<-", function(object, value) { standardGeneric("verified<-") })
#' @rdname verified-methods
#' @aliases verified<-,anlz-method
#' @export
setMethod("verified<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "verified" %in% slotNames(object) ){
              object@"verified" <- value
            } else {
              warning("verified is not in slotNames of object")
            }
            return(object)
          })
