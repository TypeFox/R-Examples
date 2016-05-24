#' @name srow_y-methods
#' @title Extract Image Attribute \code{srow_y}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{srow_y} field.  
#' @description Methods that act on the \code{srow_y} field in the
#' NIfTI/ANALYZE header.
#' @rdname srow_y-methods
#' @aliases srow_y-methods, srow_y
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
setGeneric("srow_y", function(object) standardGeneric("srow_y"))
#' @rdname srow_y-methods
#' @aliases srow_y,nifti-method
#' @export
setMethod("srow_y", "nifti", function(object) { object@"srow_y" })
#' @rdname srow_y-methods
#' @aliases srow_y<- 
#' @export
setGeneric("srow_y<-", function(object, value) { standardGeneric("srow_y<-") })
#' @rdname srow_y-methods
#' @aliases srow_y<-,nifti-method
#' @export
setMethod("srow_y<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "srow_y" %in% slotNames(object) ){
              object@"srow_y" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("srow_y <-", value))               
            } else {
              warning("srow_y is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname srow_y-methods
#' @aliases srow.y,nifti-method
#' @export
setGeneric("srow.y", function(object) standardGeneric("srow.y"))
#' @rdname srow_y-methods
#' @aliases srow.y,nifti-method
#' @export
setMethod("srow.y", "nifti", function(object) { object@"srow_y" })
#' @rdname srow_y-methods
#' @aliases srow.y<- 
#' @export
setGeneric("srow.y<-", function(object, value) { standardGeneric("srow.y<-") })
#' @rdname srow_y-methods
#' @aliases srow.y<-,nifti-method
#' @export
setMethod("srow.y<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "srow_y" %in% slotNames(object) ){
              object@"srow_y" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("srow_y <-", value))               
            } else {
              warning("srow_y is not in slotNames of object")
            }                       
            return(object)
          })
