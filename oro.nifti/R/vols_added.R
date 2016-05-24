#' @name vols_added-methods
#' @title Extract Image Attribute \code{vols_added}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{vols_added} field.  
#' @description Methods that act on the \code{vols_added} field in the
#' NIfTI/ANALYZE header.
#' @rdname vols_added-methods
#' @aliases vols_added-methods, vols_added
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
setGeneric("vols_added", function(object) standardGeneric("vols_added"))
#' @rdname vols_added-methods
#' @aliases vols_added,anlz-method
#' @export
setMethod("vols_added", "anlz", function(object) { object@"vols_added" })
#' @rdname vols_added-methods
#' @aliases vols_added<- 
#' @export
setGeneric("vols_added<-", function(object, value) { standardGeneric("vols_added<-") })
#' @rdname vols_added-methods
#' @aliases vols_added<-,anlz-method
#' @export
setMethod("vols_added<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "vols_added" %in% slotNames(object) ){
              object@"vols_added" <- value
            } else {
              warning("vols_added is not in slotNames of object")
            }
            return(object)
          })
#' @rdname vols_added-methods
#' @aliases vols.added,nifti-method
#' @export
setGeneric("vols.added", function(object) standardGeneric("vols.added"))
#' @rdname vols_added-methods
#' @aliases vols.added,anlz-method
#' @export
setMethod("vols.added", "anlz", function(object) { object@"vols_added" })
#' @rdname vols_added-methods
#' @aliases vols.added<- 
#' @export
setGeneric("vols.added<-", function(object, value) { standardGeneric("vols.added<-") })
#' @rdname vols_added-methods
#' @aliases vols.added<-,anlz-method
#' @export
setMethod("vols.added<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "vols_added" %in% slotNames(object) ){
              object@"vols_added" <- value
            } else {
              warning("vols_added is not in slotNames of object")
            }
            return(object)
          })
