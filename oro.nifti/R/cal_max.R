#' @name cal_max-methods
#' @title Extract Image Attribute \code{cal_max}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{cal_max} field.  
#' @description Methods that act on the \code{cal_max} field in the
#' NIfTI/ANALYZE header.
#' @rdname cal_max-methods
#' @aliases cal_max-methods, cal_max
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
#' cal.max(mniLR)
#' @export
setGeneric("cal_max", function(object) standardGeneric("cal_max"))
#' @rdname cal_max-methods
#' @aliases cal_max,nifti-method
#' @export
setMethod("cal_max", "nifti", function(object) { object@"cal_max" })
#' @rdname cal_max-methods
#' @aliases cal_max,anlz-method
#' @export
setMethod("cal_max", "anlz", function(object) { object@"cal_max" })
#' @rdname cal_max-methods
#' @aliases cal_max<- 
#' @export
setGeneric("cal_max<-", function(object, value) { standardGeneric("cal_max<-") })
#' @rdname cal_max-methods
#' @aliases cal_max<-,nifti-method
#' @export
setMethod("cal_max<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "cal_max" %in% slotNames(object) ){
              object@"cal_max" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("cal_max <-", value))               
            } else {
              warning("cal_max is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname cal_max-methods
#' @aliases cal_max<-,anlz-method
#' @export
setMethod("cal_max<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "cal_max" %in% slotNames(object) ){
              object@"cal_max" <- value
            } else {
              warning("cal_max is not in slotNames of object")
            }
            return(object)
          })
#' @rdname cal_max-methods
#' @aliases cal.max,nifti-method
#' @export
setGeneric("cal.max", function(object) standardGeneric("cal.max"))
#' @rdname cal_max-methods
#' @aliases cal.max,nifti-method
#' @export
setMethod("cal.max", "nifti", function(object) { object@"cal_max" })
#' @rdname cal_max-methods
#' @aliases cal.max,anlz-method
#' @export
setMethod("cal.max", "anlz", function(object) { object@"cal_max" })
#' @rdname cal_max-methods
#' @aliases cal.max<- 
#' @export
setGeneric("cal.max<-", function(object, value) { standardGeneric("cal.max<-") })
#' @rdname cal_max-methods
#' @aliases cal.max<-,nifti-method
#' @export
setMethod("cal.max<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "cal_max" %in% slotNames(object) ){
              object@"cal_max" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("cal_max <-", value))               
            } else {
              warning("cal_max is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname cal_max-methods
#' @aliases cal.max<-,anlz-method
#' @export
setMethod("cal.max<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "cal_max" %in% slotNames(object) ){
              object@"cal_max" <- value
            } else {
              warning("cal_max is not in slotNames of object")
            }
            return(object)
          })
