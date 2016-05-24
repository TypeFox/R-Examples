#' @name sizeof_hdr-methods
#' @title Extract Image Attribute \code{sizeof_hdr}
#' @docType methods 
#' @param object is an object of class \code{nifti} or \code{anlz}.
#'
#' @description Methods that act on the \code{sizeof_hdr} field in the
#' NIfTI/ANALYZE header.
#' @rdname sizeof_hdr-methods
#' @aliases sizeof_hdr-methods, sizeof_hdr
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
setGeneric("sizeof_hdr", function(object) standardGeneric("sizeof_hdr"))
#' @rdname sizeof_hdr-methods
#' @aliases sizeof_hdr,nifti-method
#' @export
setMethod("sizeof_hdr", "nifti", function(object) { object@"sizeof_hdr" })
#' @rdname sizeof_hdr-methods
#' @aliases sizeof_hdr,anlz-method
#' @export
setMethod("sizeof_hdr", "anlz", function(object) { object@"sizeof_hdr" })
#' @rdname sizeof_hdr-methods
#' @aliases sizeof.hdr,nifti-method
#' @export
setGeneric("sizeof.hdr", function(object) standardGeneric("sizeof.hdr"))
#' @rdname sizeof_hdr-methods
#' @aliases sizeof.hdr,nifti-method
#' @export
setMethod("sizeof.hdr", "nifti", function(object) { object@"sizeof_hdr" })
#' @rdname sizeof_hdr-methods
#' @aliases sizeof.hdr,anlz-method
#' @export
setMethod("sizeof.hdr", "anlz", function(object) { object@"sizeof_hdr" })
