#' @title Creates the \code{onefile} Specification for NIfTI
#' @return Object of class \code{nifti}.
#' @param img is a \code{nifti}-class object.
#' @description Changes the \code{magic} and \code{vox_offset} slots to be
#' consistent with the onefile option in \code{\link{writeNIfTI}}.  As of 
#' version 0.4.0, \code{oro.nifti} did not support the \code{"ni1"} magic type 
#' for output.
#' @author John Muschelli \email{muschellij2@@gmail.com}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @export
onefile <- function(img) {
  img@magic <- "n+1"
  img@"vox_offset" <- 352
  return(img)
}
