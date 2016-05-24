#' @title Center of Gravity for Multiple Areas
#' @description Find Center of Gravity of Each Area of Image, after thresholding
#' @param img Object of class nifti
#' @param k Minimum number of voxels for a cluster.  See \code{\link{spm_bwlabel}} 
#' @param ... Arguments passed to \code{\link{spm_bwlabel}}
#' @return Matrix of 3 columns
#' @import fslr
#' @export
group_xyz = function(img, k = 1, ...){
  
  #   stopifnot(inherits(img, "nifti"))
  les_xyz = spm_bwlabel(img, 
                        binary = FALSE, 
                        retimg = TRUE, 
                        k = k,
                        ...) 
  les_levs = sort(unique(les_xyz[les_xyz!=0]))
  if (length(les_levs) <= 1){
    return(t(xyz(les_xyz)))
  } else {
    les_xyz = t(sapply(les_levs, function(x){
      xyz(les_xyz == x)
    })) 
    return(les_xyz)
  }
  return(les_xyz)
}

  
