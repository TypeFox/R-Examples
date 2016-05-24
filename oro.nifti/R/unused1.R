#' @name unused1-methods
#' @title Extract Image Attribute \code{unused1}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{unused1} field.  
#' @description Methods that act on the \code{unused1} field in the
#' NIfTI/ANALYZE header.
#' @rdname unused1-methods
#' @aliases unused1-methods, unused1
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
setGeneric("unused1", function(object) standardGeneric("unused1"))
#' @rdname unused1-methods
#' @aliases unused1,anlz-method
#' @export
setMethod("unused1", "anlz", function(object) { object@"unused1" })
#' @rdname unused1-methods
#' @aliases unused1<- 
#' @export
setGeneric("unused1<-", function(object, value) { standardGeneric("unused1<-") })
#' @rdname unused1-methods
#' @aliases unused1<-,anlz-method
#' @export
setMethod("unused1<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "unused1" %in% slotNames(object) ){
              object@"unused1" <- value
            } else {
              warning("unused1 is not in slotNames of object")
            }
            return(object)
          })
