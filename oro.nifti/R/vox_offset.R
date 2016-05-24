#' @name vox_offset-methods
#' @title Extract Image Attribute \code{vox_offset}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{vox_offset} field.  
#' @description Methods that act on the \code{vox_offset} field in the
#' NIfTI/ANALYZE header.
#' @rdname vox_offset-methods
#' @aliases vox_offset-methods, vox_offset
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
setGeneric("vox_offset", function(object) standardGeneric("vox_offset"))
#' @rdname vox_offset-methods
#' @aliases vox_offset,nifti-method
#' @export
setMethod("vox_offset", "nifti", function(object) { object@"vox_offset" })
#' @rdname vox_offset-methods
#' @aliases vox_offset,anlz-method
#' @export
setMethod("vox_offset", "anlz", function(object) { object@"vox_offset" })
#' @rdname vox_offset-methods
#' @aliases vox_offset<- 
#' @export
setGeneric("vox_offset<-", function(object, value) { standardGeneric("vox_offset<-") })
#' @rdname vox_offset-methods
#' @aliases vox_offset<-,nifti-method
#' @export
setMethod("vox_offset<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "vox_offset" %in% slotNames(object) ){
              object@"vox_offset" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("vox_offset <-", value))               
            } else {
              warning("vox_offset is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname vox_offset-methods
#' @aliases vox_offset<-,anlz-method
#' @export
setMethod("vox_offset<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "vox_offset" %in% slotNames(object) ){
              object@"vox_offset" <- value
            } else {
              warning("vox_offset is not in slotNames of object")
            }
            return(object)
          })
#' @rdname vox_offset-methods
#' @aliases vox.offset,nifti-method
#' @export
setGeneric("vox.offset", function(object) standardGeneric("vox.offset"))
#' @rdname vox_offset-methods
#' @aliases vox.offset,nifti-method
#' @export
setMethod("vox.offset", "nifti", function(object) { object@"vox_offset" })
#' @rdname vox_offset-methods
#' @aliases vox.offset,anlz-method
#' @export
setMethod("vox.offset", "anlz", function(object) { object@"vox_offset" })
#' @rdname vox_offset-methods
#' @aliases vox.offset<- 
#' @export
setGeneric("vox.offset<-", function(object, value) { standardGeneric("vox.offset<-") })
#' @rdname vox_offset-methods
#' @aliases vox.offset<-,nifti-method
#' @export
setMethod("vox.offset<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "vox_offset" %in% slotNames(object) ){
              object@"vox_offset" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("vox_offset <-", value))               
            } else {
              warning("vox_offset is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname vox_offset-methods
#' @aliases vox.offset<-,anlz-method
#' @export
setMethod("vox.offset<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "vox_offset" %in% slotNames(object) ){
              object@"vox_offset" <- value
            } else {
              warning("vox_offset is not in slotNames of object")
            }
            return(object)
          })
