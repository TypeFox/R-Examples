#' @name datatype-methods
#' @title Extract Image Attribute \code{datatype}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{datatype} field.  
#' @description Methods that act on the \code{datatype} field in the
#' NIfTI/ANALYZE header.
#' @rdname datatype-methods
#' @aliases datatype-methods, datatype
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
setGeneric("datatype", function(object) standardGeneric("datatype"))
#' @rdname datatype-methods
#' @aliases datatype,nifti-method
#' @export
setMethod("datatype", "nifti", function(object) { object@"datatype" })
#' @rdname datatype-methods
#' @aliases datatype,anlz-method
#' @export
setMethod("datatype", "anlz", function(object) { object@"datatype" })
#' @rdname datatype-methods
#' @aliases datatype<- 
#' @export
setGeneric("datatype<-", function(object, value) { standardGeneric("datatype<-") })
#' @rdname datatype-methods
#' @aliases datatype<-,nifti-method
#' @export
setMethod("datatype<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "datatype" %in% slotNames(object) ){
              object@"datatype" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("datatype <-", value))               
            } else {
              warning("datatype is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname datatype-methods
#' @aliases datatype<-,anlz-method
#' @export
setMethod("datatype<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "datatype" %in% slotNames(object) ){
              object@"datatype" <- value
            } else {
              warning("datatype is not in slotNames of object")
            }
            return(object)
          })
