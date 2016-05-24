#' @name xyzt_units-methods
#' @title Extract Image Attribute \code{xyzt_units}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{xyzt_units} field.  
#' @description Methods that act on the \code{xyzt_units} field in the
#' NIfTI/ANALYZE header.
#' @rdname xyzt_units-methods
#' @aliases xyzt_units-methods, xyzt_units
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
setGeneric("xyzt_units", function(object) standardGeneric("xyzt_units"))
#' @rdname xyzt_units-methods
#' @aliases xyzt_units,nifti-method
#' @export
setMethod("xyzt_units", "nifti", function(object) { object@"xyzt_units" })
#' @rdname xyzt_units-methods
#' @aliases xyzt_units<- 
#' @export
setGeneric("xyzt_units<-", function(object, value) { standardGeneric("xyzt_units<-") })
#' @rdname xyzt_units-methods
#' @aliases xyzt_units<-,nifti-method
#' @export
setMethod("xyzt_units<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "xyzt_units" %in% slotNames(object) ){
              object@"xyzt_units" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("xyzt_units <-", value))               
            } else {
              warning("xyzt_units is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname xyzt_units-methods
#' @aliases xyzt.units,nifti-method
#' @export
setGeneric("xyzt.units", function(object) standardGeneric("xyzt.units"))
#' @rdname xyzt_units-methods
#' @aliases xyzt.units,nifti-method
#' @export
setMethod("xyzt.units", "nifti", function(object) { object@"xyzt_units" })
#' @rdname xyzt_units-methods
#' @aliases xyzt.units<- 
#' @export
setGeneric("xyzt.units<-", function(object, value) { standardGeneric("xyzt.units<-") })
#' @rdname xyzt_units-methods
#' @aliases xyzt.units<-,nifti-method
#' @export
setMethod("xyzt.units<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "xyzt_units" %in% slotNames(object) ){
              object@"xyzt_units" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("xyzt_units <-", value))               
            } else {
              warning("xyzt_units is not in slotNames of object")
            }                       
            return(object)
          })
