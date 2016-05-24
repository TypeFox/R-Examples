#' @name exp_time-methods
#' @title Extract Image Attribute \code{exp_time}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{exp_time} field.  
#' @description Methods that act on the \code{exp_time} field in the
#' NIfTI/ANALYZE header.
#' @rdname exp_time-methods
#' @aliases exp_time-methods, exp_time
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
setGeneric("exp_time", function(object) standardGeneric("exp_time"))
#' @rdname exp_time-methods
#' @aliases exp_time,anlz-method
#' @export
setMethod("exp_time", "anlz", function(object) { object@"exp_time" })
#' @rdname exp_time-methods
#' @aliases exp_time<- 
#' @export
setGeneric("exp_time<-", function(object, value) { standardGeneric("exp_time<-") })
#' @rdname exp_time-methods
#' @aliases exp_time<-,anlz-method
#' @export
setMethod("exp_time<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "exp_time" %in% slotNames(object) ){
              object@"exp_time" <- value
            } else {
              warning("exp_time is not in slotNames of object")
            }
            return(object)
          })
