#' @name data_type-methods
#' @title Extract Image Attribute \code{data_type}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{data_type} field.  
#' @description Methods that act on the \code{data_type} field in the
#' NIfTI/ANALYZE header.
#' @rdname data_type-methods
#' @aliases data_type-methods, data_type
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
setGeneric("data_type", function(object) standardGeneric("data_type"))
#' @rdname data_type-methods
#' @aliases data_type,nifti-method
#' @export
setMethod("data_type", "nifti", function(object) { object@"data_type" })
#' @rdname data_type-methods
#' @aliases data_type,anlz-method
#' @export
setMethod("data_type", "anlz", function(object) { object@"data_type" })
#' @rdname data_type-methods
#' @aliases data_type<- 
#' @export
setGeneric("data_type<-", function(object, value) { standardGeneric("data_type<-") })
#' @rdname data_type-methods
#' @aliases data_type<-,nifti-method
#' @export
setMethod("data_type<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "data_type" %in% slotNames(object) ){
              object@"data_type" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("data_type <-", value))               
            } else {
              warning("data_type is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname data_type-methods
#' @aliases data_type<-,anlz-method
#' @export
setMethod("data_type<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "data_type" %in% slotNames(object) ){
              object@"data_type" <- value
            } else {
              warning("data_type is not in slotNames of object")
            }
            return(object)
          })
#' @rdname data_type-methods
#' @aliases data.type,nifti-method
#' @export
setGeneric("data.type", function(object) standardGeneric("data.type"))
#' @rdname data_type-methods
#' @aliases data.type,nifti-method
#' @export
setMethod("data.type", "nifti", function(object) { object@"data_type" })
#' @rdname data_type-methods
#' @aliases data.type,anlz-method
#' @export
setMethod("data.type", "anlz", function(object) { object@"data_type" })
#' @rdname data_type-methods
#' @aliases data.type<- 
#' @export
setGeneric("data.type<-", function(object, value) { standardGeneric("data.type<-") })
#' @rdname data_type-methods
#' @aliases data.type<-,nifti-method
#' @export
setMethod("data.type<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "data_type" %in% slotNames(object) ){
              object@"data_type" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("data_type <-", value))               
            } else {
              warning("data_type is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname data_type-methods
#' @aliases data.type<-,anlz-method
#' @export
setMethod("data.type<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "data_type" %in% slotNames(object) ){
              object@"data_type" <- value
            } else {
              warning("data_type is not in slotNames of object")
            }
            return(object)
          })
