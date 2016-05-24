#' @name descrip-methods
#' @title Extract Image Attribute \code{descrip}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{descrip} field.  
#' @description Methods that act on the \code{descrip} field in the
#' NIfTI/ANALYZE header.
#' @rdname descrip-methods
#' @aliases descrip-methods, descrip
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
#' descrip(mniLR)
#' \dontrun{
#' descrip(mniLR) <- paste(descrip(mniLR), version$version.string, sep="; ")
#' descrip(mniLR)
#' }
#' @export
setGeneric("descrip", function(object) standardGeneric("descrip"))
#' @rdname descrip-methods
#' @aliases descrip,nifti-method
#' @export
setMethod("descrip", "nifti", function(object) { object@"descrip" })
#' @rdname descrip-methods
#' @aliases descrip,anlz-method
#' @export
setMethod("descrip", "anlz", function(object) { object@"descrip" })
#' @rdname descrip-methods
#' @aliases descrip<- 
#' @export
setGeneric("descrip<-", function(object, value) { standardGeneric("descrip<-") })
#' @rdname descrip-methods
#' @aliases descrip<-,nifti-method
#' @export
setMethod("descrip<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "descrip" %in% slotNames(object) ){
              object@"descrip" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("descrip <-", value))               
            } else {
              warning("descrip is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname descrip-methods
#' @aliases descrip<-,anlz-method
#' @export
setMethod("descrip<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "descrip" %in% slotNames(object) ){
              object@"descrip" <- value
            } else {
              warning("descrip is not in slotNames of object")
            }
            return(object)
          })
