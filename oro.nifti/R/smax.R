#' @name smax-methods
#' @title Extract Image Attribute \code{smax}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{smax} field.  
#' @description Methods that act on the \code{smax} field in the
#' NIfTI/ANALYZE header.
#' @rdname smax-methods
#' @aliases smax-methods, smax
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
setGeneric("smax", function(object) standardGeneric("smax"))
#' @rdname smax-methods
#' @aliases smax,anlz-method
#' @export
setMethod("smax", "anlz", function(object) { object@"smax" })
#' @rdname smax-methods
#' @aliases smax<- 
#' @export
setGeneric("smax<-", function(object, value) { standardGeneric("smax<-") })
#' @rdname smax-methods
#' @aliases smax<-,anlz-method
#' @export
setMethod("smax<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "smax" %in% slotNames(object) ){
              object@"smax" <- value
            } else {
              warning("smax is not in slotNames of object")
            }
            return(object)
          })
