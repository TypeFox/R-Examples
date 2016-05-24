
#' @title Gets Voxel Dimensions
#' @return Vector of length 3
#' @param img nifti object
#' @description Grabs the pixdim and takes the correct elements
#' @export
voxdim = function(img){
  pixdim(img)[2:4]
}

