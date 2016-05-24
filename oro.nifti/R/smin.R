#' @name smin-methods
#' @title Extract Image Attribute \code{smin}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{smin} field.  
#' @description Methods that act on the \code{smin} field in the
#' NIfTI/ANALYZE header.
#' @rdname smin-methods
#' @aliases smin-methods, smin
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
setGeneric("smin", function(object) standardGeneric("smin"))
#' @rdname smin-methods
#' @aliases smin,anlz-method
#' @export
setMethod("smin", "anlz", function(object) { object@"smin" })
#' @rdname smin-methods
#' @aliases smin<- 
#' @export
setGeneric("smin<-", function(object, value) { standardGeneric("smin<-") })
#' @rdname smin-methods
#' @aliases smin<-,anlz-method
#' @export
setMethod("smin<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "smin" %in% slotNames(object) ){
              object@"smin" <- value
            } else {
              warning("smin is not in slotNames of object")
            }
            return(object)
          })
