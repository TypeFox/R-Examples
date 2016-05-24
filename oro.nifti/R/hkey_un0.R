#' @name hkey_un0-methods
#' @title Extract Image Attribute \code{hkey_un0}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{hkey_un0} field.  
#' @description Methods that act on the \code{hkey_un0} field in the
#' NIfTI/ANALYZE header.
#' @rdname hkey_un0-methods
#' @aliases hkey_un0-methods, hkey_un0
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
setGeneric("hkey_un0", function(object) standardGeneric("hkey_un0"))
#' @rdname hkey_un0-methods
#' @aliases hkey_un0,anlz-method
#' @export
setMethod("hkey_un0", "anlz", function(object) { object@"hkey_un0" })
#' @rdname hkey_un0-methods
#' @aliases hkey_un0<- 
#' @export
setGeneric("hkey_un0<-", function(object, value) { standardGeneric("hkey_un0<-") })
#' @rdname hkey_un0-methods
#' @aliases hkey_un0<-,anlz-method
#' @export
setMethod("hkey_un0<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "hkey_un0" %in% slotNames(object) ){
              object@"hkey_un0" <- value
            } else {
              warning("hkey_un0 is not in slotNames of object")
            }
            return(object)
          })
#' @rdname hkey_un0-methods
#' @aliases hkey.un0,nifti-method
#' @export
setGeneric("hkey.un0", function(object) standardGeneric("hkey.un0"))
#' @rdname hkey_un0-methods
#' @aliases hkey.un0,anlz-method
#' @export
setMethod("hkey.un0", "anlz", function(object) { object@"hkey_un0" })
#' @rdname hkey_un0-methods
#' @aliases hkey.un0<- 
#' @export
setGeneric("hkey.un0<-", function(object, value) { standardGeneric("hkey.un0<-") })
#' @rdname hkey_un0-methods
#' @aliases hkey.un0<-,anlz-method
#' @export
setMethod("hkey.un0<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "hkey_un0" %in% slotNames(object) ){
              object@"hkey_un0" <- value
            } else {
              warning("hkey_un0 is not in slotNames of object")
            }
            return(object)
          })
