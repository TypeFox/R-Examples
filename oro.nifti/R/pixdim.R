#' @name pixdim-methods
#' @title Extract Image Attribute \code{pixdim}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{pixdim} field.  
#' @description Methods that act on the \code{pixdim} field in the
#' NIfTI/ANALYZE header.
#' @rdname pixdim-methods
#' @aliases pixdim-methods, pixdim
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
#' "mniLR.nii.gz")
#' download.file(url, urlfile, quiet=TRUE)
#' }
#' urlfile <- file.path(system.file("nifti", package="oro.nifti"),
#'                      "mniLR.nii.gz")
#' mniLR <- readNIfTI(urlfile)
#' pixdim(mniLR)
#' @export
setGeneric("pixdim", function(object) standardGeneric("pixdim"))
#' @rdname pixdim-methods
#' @aliases pixdim,nifti-method
#' @export
setMethod("pixdim", "nifti", function(object) { object@"pixdim" })
#' @rdname pixdim-methods
#' @aliases pixdim,anlz-method
#' @export
setMethod("pixdim", "anlz", function(object) { object@"pixdim" })
#' @rdname pixdim-methods
#' @aliases pixdim<- 
#' @export
setGeneric("pixdim<-", function(object, value) { standardGeneric("pixdim<-") })
#' @rdname pixdim-methods
#' @aliases pixdim<-,nifti-method
#' @export
setMethod("pixdim<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "pixdim" %in% slotNames(object) ){
              object@"pixdim" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("pixdim <-", value))               
            } else {
              warning("pixdim is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname pixdim-methods
#' @aliases pixdim<-,anlz-method
#' @export
setMethod("pixdim<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "pixdim" %in% slotNames(object) ){
              object@"pixdim" <- value
            } else {
              warning("pixdim is not in slotNames of object")
            }
            return(object)
          })
