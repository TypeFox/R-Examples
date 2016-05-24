#' @name orient-methods
#' @title Extract Image Attribute \code{orient}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{orient} field.  
#' @description Methods that act on the \code{orient} field in the
#' NIfTI/ANALYZE header.
#' @rdname orient-methods
#' @aliases orient-methods, orient
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
setGeneric("orient", function(object) standardGeneric("orient"))
#' @rdname orient-methods
#' @aliases orient,anlz-method
#' @export
setMethod("orient", "anlz", function(object) { object@"orient" })
#' @rdname orient-methods
#' @aliases orient<- 
#' @export
setGeneric("orient<-", function(object, value) { standardGeneric("orient<-") })
#' @rdname orient-methods
#' @aliases orient<-,anlz-method
#' @export
setMethod("orient<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "orient" %in% slotNames(object) ){
              object@"orient" <- value
            } else {
              warning("orient is not in slotNames of object")
            }
            return(object)
          })
