#' @name slice_end-methods
#' @title Extract Image Attribute \code{slice_end}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{slice_end} field.  
#' @description Methods that act on the \code{slice_end} field in the
#' NIfTI/ANALYZE header.
#' @rdname slice_end-methods
#' @aliases slice_end-methods, slice_end
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
setGeneric("slice_end", function(object) standardGeneric("slice_end"))
#' @rdname slice_end-methods
#' @aliases slice_end,nifti-method
#' @export
setMethod("slice_end", "nifti", function(object) { object@"slice_end" })
#' @rdname slice_end-methods
#' @aliases slice_end<- 
#' @export
setGeneric("slice_end<-", function(object, value) { standardGeneric("slice_end<-") })
#' @rdname slice_end-methods
#' @aliases slice_end<-,nifti-method
#' @export
setMethod("slice_end<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "slice_end" %in% slotNames(object) ){
              object@"slice_end" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("slice_end <-", value))               
            } else {
              warning("slice_end is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname slice_end-methods
#' @aliases slice.end,nifti-method
#' @export
setGeneric("slice.end", function(object) standardGeneric("slice.end"))
#' @rdname slice_end-methods
#' @aliases slice.end,nifti-method
#' @export
setMethod("slice.end", "nifti", function(object) { object@"slice_end" })
#' @rdname slice_end-methods
#' @aliases slice.end<- 
#' @export
setGeneric("slice.end<-", function(object, value) { standardGeneric("slice.end<-") })
#' @rdname slice_end-methods
#' @aliases slice.end<-,nifti-method
#' @export
setMethod("slice.end<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "slice_end" %in% slotNames(object) ){
              object@"slice_end" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("slice_end <-", value))               
            } else {
              warning("slice_end is not in slotNames of object")
            }                       
            return(object)
          })
