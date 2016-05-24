#' @name aux_file-methods
#' @title Extract Image Attribute \code{aux_file}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#' @param value is the value to assign to the \code{aux_file} field.  
#' @description Methods that act on the \code{aux_file} field in the
#' NIfTI/ANALYZE header.
#' @rdname aux_file-methods
#' @aliases aux_file-methods, aux_file
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
#' url <- "http://nifti.nimh.nih.gov/nifti-1/data/avg152T1_RL_nifti.nii.gz"
#' urlfile <- file.path(system.file("nifti", package="oro.nifti"),
#'                     "mniRL.nii.gz")
#' download.file(url, urlfile, quiet=TRUE)
#' }
#' options("niftiAuditTrail"=FALSE)
#' 
#' urlfile <- file.path(system.file("nifti", package="oro.nifti"),
#'                      "mniRL.nii.gz")
#' mniRL <- readNIfTI(urlfile)
#' aux.file(mniRL)
#' aux.file(mniRL) <- "avg152T1_RL_nifti"
#' aux.file(mniRL)
#' @export
setGeneric("aux_file", function(object) standardGeneric("aux_file"))
#' @rdname aux_file-methods
#' @aliases aux_file,nifti-method
#' @export
setMethod("aux_file", "nifti", function(object) { object@"aux_file" })
#' @rdname aux_file-methods
#' @aliases aux_file,anlz-method
#' @export
setMethod("aux_file", "anlz", function(object) { object@"aux_file" })
#' @rdname aux_file-methods
#' @aliases aux_file<- 
#' @export
setGeneric("aux_file<-", function(object, value) { standardGeneric("aux_file<-") })
#' @rdname aux_file-methods
#' @aliases aux_file<-,nifti-method
#' @export
setMethod("aux_file<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "aux_file" %in% slotNames(object) ){
              object@"aux_file" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("aux_file <-", value))               
            } else {
              warning("aux_file is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname aux_file-methods
#' @aliases aux_file<-,anlz-method
#' @export
setMethod("aux_file<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "aux_file" %in% slotNames(object) ){
              object@"aux_file" <- value
            } else {
              warning("aux_file is not in slotNames of object")
            }
            return(object)
          })
#' @rdname aux_file-methods
#' @aliases aux.file,nifti-method
#' @export
setGeneric("aux.file", function(object) standardGeneric("aux.file"))
#' @rdname aux_file-methods
#' @aliases aux.file,nifti-method
#' @export
setMethod("aux.file", "nifti", function(object) { object@"aux_file" })
#' @rdname aux_file-methods
#' @aliases aux.file,anlz-method
#' @export
setMethod("aux.file", "anlz", function(object) { object@"aux_file" })
#' @rdname aux_file-methods
#' @aliases aux.file<- 
#' @export
setGeneric("aux.file<-", function(object, value) { standardGeneric("aux.file<-") })
#' @rdname aux_file-methods
#' @aliases aux.file<-,nifti-method
#' @export
setMethod("aux.file<-", 
          signature(object="nifti"), 
          function(object, value) { 
            if ( "aux_file" %in% slotNames(object) ){
              object@"aux_file" <- value
              audit.trail(object) <-
                niftiAuditTrailEvent(object, "modification", match.call(),
                                     paste("aux_file <-", value))               
            } else {
              warning("aux_file is not in slotNames of object")
            }                       
            return(object)
          })
#' @rdname aux_file-methods
#' @aliases aux.file<-,anlz-method
#' @export
setMethod("aux.file<-", 
          signature(object="anlz"), 
          function(object, value) { 
            if ( "aux_file" %in% slotNames(object) ){
              object@"aux_file" <- value
            } else {
              warning("aux_file is not in slotNames of object")
            }
            return(object)
          })
