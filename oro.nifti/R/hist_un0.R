#' @name hist_un0-methods
#' @title Extract Image Attribute \code{hist_un0}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{hist_un0} field.  
#' @description Methods that act on the \code{hist_un0} field in the
#' NIfTI/ANALYZE header.
#' @rdname hist_un0-methods
#' @aliases hist_un0-methods, hist_un0
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
setGeneric("hist_un0", function(object) standardGeneric("hist_un0"))
#' @rdname hist_un0-methods
#' @aliases hist_un0,anlz-method
#' @export
setMethod("hist_un0", "anlz", function(object) { object@"hist_un0" })
#' @rdname hist_un0-methods
#' @aliases hist_un0<- 
#' @export
setGeneric("hist_un0<-", function(object, value) { standardGeneric("hist_un0<-") })
#' @rdname hist_un0-methods
#' @aliases hist_un0<-,anlz-method
#' @export
setMethod("hist_un0<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "hist_un0" %in% slotNames(object) ){
              object@"hist_un0" <- value
            } else {
              warning("hist_un0 is not in slotNames of object")
            }
            return(object)
          })
