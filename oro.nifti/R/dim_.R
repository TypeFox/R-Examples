#' @name dim_-methods
#' @title Extract Image Attribute \code{dim_}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{dim_} field.  
#' @description Methods that act on the \code{dim_} field in the
#' NIfTI/ANALYZE header.
#' @rdname dim_-methods
#' @aliases dim_-methods, dim_
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
setGeneric("dim_", function(object) standardGeneric("dim_"))
#' @rdname dim_-methods
#' @aliases dim_,nifti-method
#' @export
setMethod("dim_", "nifti", function(object) { object@"dim_" })
#' @rdname dim_-methods
#' @aliases dim_,anlz-method
#' @export
setMethod("dim_", "anlz", function(object) { object@"dim_" })
#' @rdname dim_-methods
#' @aliases dim_<- 
#' @export
setGeneric("dim_<-", function(object, value) { standardGeneric("dim_<-") })
#' @rdname dim_-methods
#' @aliases dim_<-,nifti-method
#' @export
setMethod("dim_<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "dim_" %in% slotNames(object) ){
              object@"dim_" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("dim_ <-", value))               
            } else {
              warning("dim_ is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname dim_-methods
#' @aliases dim_<-,anlz-method
#' @export
setMethod("dim_<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "dim_" %in% slotNames(object) ){
              object@"dim_" <- value
            } else {
              warning("dim_ is not in slotNames of object")
            }
            return(object)
          })
