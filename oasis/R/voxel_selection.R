#' @title Voxel Selection Procedure
#' @description This function creates a binary mask for the voxel selection 
#' procedure for OASIS.
#' @param flair FLAIR volume of class \code{\link{nifti}}
#' @param brain_mask brain mask of class \code{\link{nifti}}
#' @param cutoff the percentile cutoff for the thresholding, 
#' passed to \code{\link{quantile}}
#' @return Returns the voxel selection mask as an object of class \code{\link{nifti}}.  
#' @examples
#' library(oasis)
#' niis = tempfile(fileext = ".nii.gz")
#' if (require(downloader)) {
#'  download(paste0("https://dl.dropbox.com/u/2785709/brainder", 
#'          "/software/flair/templates/GG-853-FLAIR-2.0mm.nii.gz"), 
#'          destfile = niis)
#'  
#'  flair <- readnii(niis) 
#'  if (have.fsl()) {
#'    brain_mask = fslbet(niis) > 0
#'  } else {
#'    ind = list(c(10L, 81L), c(12L, 101L), c(3L, 78L))
#'    all.ind = lapply(ind, function(x) seq(x[1], x[2]))
#'    brain_mask = niftiarr(flair, 0)
#'    eg = expand.grid(all.ind)
#'    eg = as.matrix(eg)
#'    brain_mask[eg] = 1
#'  }
#'  voxel_selection_mask <- voxel_selection(flair, 
#'    brain_mask, cutoff = .85) 
#' }
#' @import oro.nifti
#' @export
voxel_selection <- function(flair, ##the flair volume
                            brain_mask, ## a brain mask for the flair 
                            cutoff ## the percentile cutoff for the thresholding 
                            ){
  stopifnot(length(cutoff) == 1)
  arr = array(brain_mask, dim = dim(brain_mask))
  class(arr) = "numeric"
  stopifnot(all(c(arr) %in% c(0, 1)))
  ##find the value to threshold the flair volume at 
  cutpoint <- quantile(flair[brain_mask == 1], cutoff)
  outmask <-  flair > cutpoint
  outmask <- as.nifti(outmask * brain_mask) 
  outmask <- datatyper(outmask, trybyte = TRUE)
  ##return the binary mask of the flair values above the cutpoint 
  return(outmask) 
}