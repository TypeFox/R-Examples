#' @name dim_info-methods
#' @title Extract Image Attribute \code{dim_info}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{dim_info} field.  
#' @description Methods that act on the \code{dim_info} field in the
#' NIfTI/ANALYZE header.
#' @rdname dim_info-methods
#' @aliases dim_info-methods, dim_info
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
setGeneric("dim_info", function(object) standardGeneric("dim_info"))
#' @rdname dim_info-methods
#' @aliases dim_info,nifti-method
#' @export
setMethod("dim_info", "nifti", function(object) { object@"dim_info" })
#' @rdname dim_info-methods
#' @aliases dim_info<- 
#' @export
setGeneric("dim_info<-", function(object, value) { standardGeneric("dim_info<-") })
#' @rdname dim_info-methods
#' @aliases dim_info<-,nifti-method
#' @export
setMethod("dim_info<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "dim_info" %in% slotNames(object) ){
              object@"dim_info" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("dim_info <-", value))               
            } else {
              warning("dim_info is not in slotNames of object")
            }                       
            return(object)
          })
