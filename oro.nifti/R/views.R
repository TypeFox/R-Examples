#' @name views-methods
#' @title Extract Image Attribute \code{views}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{views} field.  
#' @description Methods that act on the \code{views} field in the
#' NIfTI/ANALYZE header.
#' @rdname views-methods
#' @aliases views-methods, views
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
setGeneric("views", function(object) standardGeneric("views"))
#' @rdname views-methods
#' @aliases views,anlz-method
#' @export
setMethod("views", "anlz", function(object) { object@"views" })
#' @rdname views-methods
#' @aliases views<- 
#' @export
setGeneric("views<-", function(object, value) { standardGeneric("views<-") })
#' @rdname views-methods
#' @aliases views<-,anlz-method
#' @export
setMethod("views<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "views" %in% slotNames(object) ){
              object@"views" <- value
            } else {
              warning("views is not in slotNames of object")
            }
            return(object)
          })
