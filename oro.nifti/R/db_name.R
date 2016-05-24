#' @name db_name-methods
#' @title Extract Image Attribute \code{db_name}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{db_name} field.  
#' @description Methods that act on the \code{db_name} field in the
#' NIfTI/ANALYZE header.
#' @rdname db_name-methods
#' @aliases db_name-methods, db_name
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
setGeneric("db_name", function(object) standardGeneric("db_name"))
#' @rdname db_name-methods
#' @aliases db_name,nifti-method
#' @export
setMethod("db_name", "nifti", function(object) { object@"db_name" })
#' @rdname db_name-methods
#' @aliases db_name,anlz-method
#' @export
setMethod("db_name", "anlz", function(object) { object@"db_name" })
#' @rdname db_name-methods
#' @aliases db_name<- 
#' @export
setGeneric("db_name<-", function(object, value) { standardGeneric("db_name<-") })
#' @rdname db_name-methods
#' @aliases db_name<-,nifti-method
#' @export
setMethod("db_name<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "db_name" %in% slotNames(object) ){
              object@"db_name" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("db_name <-", value))               
            } else {
              warning("db_name is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname db_name-methods
#' @aliases db_name<-,anlz-method
#' @export
setMethod("db_name<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "db_name" %in% slotNames(object) ){
              object@"db_name" <- value
            } else {
              warning("db_name is not in slotNames of object")
            }
            return(object)
          })
#' @rdname db_name-methods
#' @aliases db.name,nifti-method
#' @export
setGeneric("db.name", function(object) standardGeneric("db.name"))
#' @rdname db_name-methods
#' @aliases db.name,nifti-method
#' @export
setMethod("db.name", "nifti", function(object) { object@"db_name" })
#' @rdname db_name-methods
#' @aliases db.name,anlz-method
#' @export
setMethod("db.name", "anlz", function(object) { object@"db_name" })
#' @rdname db_name-methods
#' @aliases db.name<- 
#' @export
setGeneric("db.name<-", function(object, value) { standardGeneric("db.name<-") })
#' @rdname db_name-methods
#' @aliases db.name<-,nifti-method
#' @export
setMethod("db.name<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "db_name" %in% slotNames(object) ){
              object@"db_name" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("db_name <-", value))               
            } else {
              warning("db_name is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname db_name-methods
#' @aliases db.name<-,anlz-method
#' @export
setMethod("db.name<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "db_name" %in% slotNames(object) ){
              object@"db_name" <- value
            } else {
              warning("db_name is not in slotNames of object")
            }
            return(object)
          })
