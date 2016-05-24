#' @name toffset-methods
#' @title Extract Image Attribute \code{toffset}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{toffset} field.  
#' @description Methods that act on the \code{toffset} field in the
#' NIfTI/ANALYZE header.
#' @rdname toffset-methods
#' @aliases toffset-methods, toffset
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
setGeneric("toffset", function(object) standardGeneric("toffset"))
#' @rdname toffset-methods
#' @aliases toffset,nifti-method
#' @export
setMethod("toffset", "nifti", function(object) { object@"toffset" })
#' @rdname toffset-methods
#' @aliases toffset<- 
#' @export
setGeneric("toffset<-", function(object, value) { standardGeneric("toffset<-") })
#' @rdname toffset-methods
#' @aliases toffset<-,nifti-method
#' @export
setMethod("toffset<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "toffset" %in% slotNames(object) ){
              object@"toffset" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("toffset <-", value))               
            } else {
              warning("toffset is not in slotNames of object")
            }                       
            return(object)
          })
