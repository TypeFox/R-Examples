#' @name cal_units-methods
#' @title Extract Image Attribute \code{cal_units}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{cal_units} field.  
#' @description Methods that act on the \code{cal_units} field in the
#' NIfTI/ANALYZE header.
#' @rdname cal_units-methods
#' @aliases cal_units-methods, cal_units
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
setGeneric("cal_units", function(object) standardGeneric("cal_units"))
#' @rdname cal_units-methods
#' @aliases cal_units,anlz-method
#' @export
setMethod("cal_units", "anlz", function(object) { object@"cal_units" })
#' @rdname cal_units-methods
#' @aliases cal_units<- 
#' @export
setGeneric("cal_units<-", function(object, value) { standardGeneric("cal_units<-") })
#' @rdname cal_units-methods
#' @aliases cal_units<-,anlz-method
#' @export
setMethod("cal_units<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "cal_units" %in% slotNames(object) ){
              object@"cal_units" <- value
            } else {
              warning("cal_units is not in slotNames of object")
            }
            return(object)
          })
#' @rdname cal_units-methods
#' @aliases cal.units,nifti-method
#' @export
setGeneric("cal.units", function(object) standardGeneric("cal.units"))
#' @rdname cal_units-methods
#' @aliases cal.units,anlz-method
#' @export
setMethod("cal.units", "anlz", function(object) { object@"cal_units" })
#' @rdname cal_units-methods
#' @aliases cal.units<- 
#' @export
setGeneric("cal.units<-", function(object, value) { standardGeneric("cal.units<-") })
#' @rdname cal_units-methods
#' @aliases cal.units<-,anlz-method
#' @export
setMethod("cal.units<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "cal_units" %in% slotNames(object) ){
              object@"cal_units" <- value
            } else {
              warning("cal_units is not in slotNames of object")
            }
            return(object)
          })
