#' @name generated-methods
#' @title Extract Image Attribute \code{generated}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{generated} field.  
#' @description Methods that act on the \code{generated} field in the
#' NIfTI/ANALYZE header.
#' @rdname generated-methods
#' @aliases generated-methods, generated
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
setGeneric("generated", function(object) standardGeneric("generated"))
#' @rdname generated-methods
#' @aliases generated,anlz-method
#' @export
setMethod("generated", "anlz", function(object) { object@"generated" })
#' @rdname generated-methods
#' @aliases generated<- 
#' @export
setGeneric("generated<-", function(object, value) { standardGeneric("generated<-") })
#' @rdname generated-methods
#' @aliases generated<-,anlz-method
#' @export
setMethod("generated<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "generated" %in% slotNames(object) ){
              object@"generated" <- value
            } else {
              warning("generated is not in slotNames of object")
            }
            return(object)
          })
