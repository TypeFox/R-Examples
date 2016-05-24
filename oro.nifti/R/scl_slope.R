#' @name scl_slope-methods
#' @title Extract Image Attribute \code{scl_slope}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{scl_slope} field.  
#' @description Methods that act on the \code{scl_slope} field in the
#' NIfTI/ANALYZE header.
#' @rdname scl_slope-methods
#' @aliases scl_slope-methods, scl_slope
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
setGeneric("scl_slope", function(object) standardGeneric("scl_slope"))
#' @rdname scl_slope-methods
#' @aliases scl_slope,nifti-method
#' @export
setMethod("scl_slope", "nifti", function(object) { object@"scl_slope" })
#' @rdname scl_slope-methods
#' @aliases scl_slope<- 
#' @export
setGeneric("scl_slope<-", function(object, value) { standardGeneric("scl_slope<-") })
#' @rdname scl_slope-methods
#' @aliases scl_slope<-,nifti-method
#' @export
setMethod("scl_slope<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "scl_slope" %in% slotNames(object) ){
              object@"scl_slope" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("scl_slope <-", value))               
            } else {
              warning("scl_slope is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname scl_slope-methods
#' @aliases scl.slope,nifti-method
#' @export
setGeneric("scl.slope", function(object) standardGeneric("scl.slope"))
#' @rdname scl_slope-methods
#' @aliases scl.slope,nifti-method
#' @export
setMethod("scl.slope", "nifti", function(object) { object@"scl_slope" })
#' @rdname scl_slope-methods
#' @aliases scl.slope<- 
#' @export
setGeneric("scl.slope<-", function(object, value) { standardGeneric("scl.slope<-") })
#' @rdname scl_slope-methods
#' @aliases scl.slope<-,nifti-method
#' @export
setMethod("scl.slope<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "scl_slope" %in% slotNames(object) ){
              object@"scl_slope" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("scl_slope <-", value))               
            } else {
              warning("scl_slope is not in slotNames of object")
            }                       
            return(object)
          })
