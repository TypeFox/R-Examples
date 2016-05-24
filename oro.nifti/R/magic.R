#' @name magic-methods
#' @title Extract Image Attribute \code{magic}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{magic} field.  
#' @description Methods that act on the \code{magic} field in the
#' NIfTI/ANALYZE header.
#' @rdname magic-methods
#' @aliases magic-methods, magic
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
setGeneric("magic", function(object) standardGeneric("magic"))
#' @rdname magic-methods
#' @aliases magic,nifti-method
#' @export
setMethod("magic", "nifti", function(object) { object@"magic" })
#' @rdname magic-methods
#' @aliases magic<- 
#' @export
setGeneric("magic<-", function(object, value) { standardGeneric("magic<-") })
#' @rdname magic-methods
#' @aliases magic<-,nifti-method
#' @export
setMethod("magic<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "magic" %in% slotNames(object) ){
              object@"magic" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("magic <-", value))               
            } else {
              warning("magic is not in slotNames of object")
            }                       
            return(object)
          })
