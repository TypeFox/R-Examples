#' @name vox_units-methods
#' @title Extract Image Attribute \code{vox_units}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{vox_units} field.  
#' @description Methods that act on the \code{vox_units} field in the
#' NIfTI/ANALYZE header.
#' @rdname vox_units-methods
#' @aliases vox_units-methods, vox_units
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
setGeneric("vox_units", function(object) standardGeneric("vox_units"))
#' @rdname vox_units-methods
#' @aliases vox_units,anlz-method
#' @export
setMethod("vox_units", "anlz", function(object) { object@"vox_units" })
#' @rdname vox_units-methods
#' @aliases vox_units<- 
#' @export
setGeneric("vox_units<-", function(object, value) { standardGeneric("vox_units<-") })
#' @rdname vox_units-methods
#' @aliases vox_units<-,anlz-method
#' @export
setMethod("vox_units<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "vox_units" %in% slotNames(object) ){
              object@"vox_units" <- value
            } else {
              warning("vox_units is not in slotNames of object")
            }
            return(object)
          })
#' @rdname vox_units-methods
#' @aliases vox.units,nifti-method
#' @export
setGeneric("vox.units", function(object) standardGeneric("vox.units"))
#' @rdname vox_units-methods
#' @aliases vox.units,anlz-method
#' @export
setMethod("vox.units", "anlz", function(object) { object@"vox_units" })
#' @rdname vox_units-methods
#' @aliases vox.units<- 
#' @export
setGeneric("vox.units<-", function(object, value) { standardGeneric("vox.units<-") })
#' @rdname vox_units-methods
#' @aliases vox.units<-,anlz-method
#' @export
setMethod("vox.units<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "vox_units" %in% slotNames(object) ){
              object@"vox_units" <- value
            } else {
              warning("vox_units is not in slotNames of object")
            }
            return(object)
          })
