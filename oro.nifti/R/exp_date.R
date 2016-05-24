#' @name exp_date-methods
#' @title Extract Image Attribute \code{exp_date}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{exp_date} field.  
#' @description Methods that act on the \code{exp_date} field in the
#' NIfTI/ANALYZE header.
#' @rdname exp_date-methods
#' @aliases exp_date-methods, exp_date
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
setGeneric("exp_date", function(object) standardGeneric("exp_date"))
#' @rdname exp_date-methods
#' @aliases exp_date,anlz-method
#' @export
setMethod("exp_date", "anlz", function(object) { object@"exp_date" })
#' @rdname exp_date-methods
#' @aliases exp_date<- 
#' @export
setGeneric("exp_date<-", function(object, value) { standardGeneric("exp_date<-") })
#' @rdname exp_date-methods
#' @aliases exp_date<-,anlz-method
#' @export
setMethod("exp_date<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "exp_date" %in% slotNames(object) ){
              object@"exp_date" <- value
            } else {
              warning("exp_date is not in slotNames of object")
            }
            return(object)
          })
