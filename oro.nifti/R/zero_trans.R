#' @title Change Intercept to 0 and Slope to 1 in NIfTI Object
#' @return An object of the same type passed.  
#' @param img is a \code{nifti} object (or character of filename).  If an 
#' \code{anlz} object is passed, the unaltered \code{anlz} object is returned.
#' @description Forces image \code{scl_slope} to 1 and \code{scl_inter}
#' to be 0 of slots of class \code{nifti}.  This is so that when images are 
#' rendered/written, the values correspond to those in the array (stored in the
#' \code{.Data} slot) and are not scaled.
#' @author John Muschelli \email{muschellij2@@gmail.com}
#' @name resetSlopeIntercept
#' @rdname zero_trans
#' @export
resetSlopeIntercept <- function(img) {
  if (is.nifti(img)) {
    scl.slope(img) <- 1
    scl.inter(img) <- 0
  }
  return(img)
}
#' @rdname zero_trans
#' @export
zero_trans <- function(img) {
  resetSlopeIntercept(img = img)
}
