#' @name dim_un0-methods
#' @title Extract Image Attribute \code{dim_un0}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{dim_un0} field.  
#' @description Methods that act on the \code{dim_un0} field in the
#' NIfTI/ANALYZE header.
#' @rdname dim_un0-methods
#' @aliases dim_un0-methods, dim_un0
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
setGeneric("dim_un0", function(object) standardGeneric("dim_un0"))
#' @rdname dim_un0-methods
#' @aliases dim_un0,anlz-method
#' @export
setMethod("dim_un0", "anlz", function(object) { object@"dim_un0" })
#' @rdname dim_un0-methods
#' @aliases dim_un0<- 
#' @export
setGeneric("dim_un0<-", function(object, value) { standardGeneric("dim_un0<-") })
#' @rdname dim_un0-methods
#' @aliases dim_un0<-,anlz-method
#' @export
setMethod("dim_un0<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "dim_un0" %in% slotNames(object) ){
              object@"dim_un0" <- value
            } else {
              warning("dim_un0 is not in slotNames of object")
            }
            return(object)
          })
