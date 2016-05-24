#' @name cal_min-methods
#' @title Extract Image Attribute \code{cal_min}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{cal_min} field.  
#' @description Methods that act on the \code{cal_min} field in the
#' NIfTI/ANALYZE header.
#' @rdname cal_min-methods
#' @aliases cal_min-methods, cal_min
#' @details See documentation on the ANALYZE and/or NIfTI data standards for
#' more details.
#' @author John Muschelli \email{muschellij2@@gmail.com},\cr
#' Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @references
#' ANALYZE 7.5\cr
#' \url{https://rportal.mayo.edu/bir/ANALYZE75.pdf}\cr
#' NIfTI-1\cr
#' \url{http://nifti.nimh.nih.gov/}
#' @examples \dontrun{
#' url <- "http://nifti.nimh.nih.gov/nifti-1/data/avg152T1_LR_nifti.nii.gz"
#' urlfile <- file.path(system.file("nifti", package="oro.nifti"),
#'                      "mniLR.nii.gz")
#' download.file(url, urlfile, quiet=TRUE)
#' }
#' urlfile <- file.path(system.file("nifti", package="oro.nifti"),
#'                      "mniLR.nii.gz")
#' mniLR <- readNIfTI(urlfile)
#' cal.min(mniLR)
#' @export
setGeneric("cal_min", function(object) standardGeneric("cal_min"))
#' @rdname cal_min-methods
#' @aliases cal_min,nifti-method
#' @export
setMethod("cal_min", "nifti", function(object) { object@"cal_min" })
#' @rdname cal_min-methods
#' @aliases cal_min,anlz-method
#' @export
setMethod("cal_min", "anlz", function(object) { object@"cal_min" })
#' @rdname cal_min-methods
#' @aliases cal_min<- 
#' @export
setGeneric("cal_min<-", function(object, value) { standardGeneric("cal_min<-") })
#' @rdname cal_min-methods
#' @aliases cal_min<-,nifti-method
#' @export
setMethod("cal_min<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "cal_min" %in% slotNames(object) ){
              object@"cal_min" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("cal_min <-", value))               
            } else {
              warning("cal_min is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname cal_min-methods
#' @aliases cal_min<-,anlz-method
#' @export
setMethod("cal_min<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "cal_min" %in% slotNames(object) ){
              object@"cal_min" <- value
            } else {
              warning("cal_min is not in slotNames of object")
            }
            return(object)
          })
#' @rdname cal_min-methods
#' @aliases cal.min,nifti-method
#' @export
setGeneric("cal.min", function(object) standardGeneric("cal.min"))
#' @rdname cal_min-methods
#' @aliases cal.min,nifti-method
#' @export
setMethod("cal.min", "nifti", function(object) { object@"cal_min" })
#' @rdname cal_min-methods
#' @aliases cal.min,anlz-method
#' @export
setMethod("cal.min", "anlz", function(object) { object@"cal_min" })
#' @rdname cal_min-methods
#' @aliases cal.min<- 
#' @export
setGeneric("cal.min<-", function(object, value) { standardGeneric("cal.min<-") })
#' @rdname cal_min-methods
#' @aliases cal.min<-,nifti-method
#' @export
setMethod("cal.min<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "cal_min" %in% slotNames(object) ){
              object@"cal_min" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("cal_min <-", value))               
            } else {
              warning("cal_min is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname cal_min-methods
#' @aliases cal.min<-,anlz-method
#' @export
setMethod("cal.min<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "cal_min" %in% slotNames(object) ){
              object@"cal_min" <- value
            } else {
              warning("cal_min is not in slotNames of object")
            }
            return(object)
          })
