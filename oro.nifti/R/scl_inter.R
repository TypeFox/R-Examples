#' @name scl_inter-methods
#' @title Extract Image Attribute \code{scl_inter}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{scl_inter} field.  
#' @description Methods that act on the \code{scl_inter} field in the
#' NIfTI/ANALYZE header.
#' @rdname scl_inter-methods
#' @aliases scl_inter-methods, scl_inter
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
setGeneric("scl_inter", function(object) standardGeneric("scl_inter"))
#' @rdname scl_inter-methods
#' @aliases scl_inter,nifti-method
#' @export
setMethod("scl_inter", "nifti", function(object) { object@"scl_inter" })
#' @rdname scl_inter-methods
#' @aliases scl_inter<- 
#' @export
setGeneric("scl_inter<-", function(object, value) { standardGeneric("scl_inter<-") })
#' @rdname scl_inter-methods
#' @aliases scl_inter<-,nifti-method
#' @export
setMethod("scl_inter<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "scl_inter" %in% slotNames(object) ){
              object@"scl_inter" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("scl_inter <-", value))               
            } else {
              warning("scl_inter is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname scl_inter-methods
#' @aliases scl.inter,nifti-method
#' @export
setGeneric("scl.inter", function(object) standardGeneric("scl.inter"))
#' @rdname scl_inter-methods
#' @aliases scl.inter,nifti-method
#' @export
setMethod("scl.inter", "nifti", function(object) { object@"scl_inter" })
#' @rdname scl_inter-methods
#' @aliases scl.inter<- 
#' @export
setGeneric("scl.inter<-", function(object, value) { standardGeneric("scl.inter<-") })
#' @rdname scl_inter-methods
#' @aliases scl.inter<-,nifti-method
#' @export
setMethod("scl.inter<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "scl_inter" %in% slotNames(object) ){
              object@"scl_inter" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("scl_inter <-", value))               
            } else {
              warning("scl_inter is not in slotNames of object")
            }                       
            return(object)
          })
