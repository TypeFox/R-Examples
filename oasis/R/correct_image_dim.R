#' @title Image Dimension Correction
#' @description This function takes an image and drops dimensions 
#' until the volume is a user specified dimension. 
#' @param image volume of class \code{\link{nifti}}
#' @param dim scalar value of desired image dimension
#' @import fslr
#' @return Returns a volume of class \code{\link{nifti}} of desired dimension.  
#' @examples \dontrun{
#' library(fslr)
#' flair <- readnii('path/to/flair', reorient = FALSE) 
#' flair <- correct_image_dim(flair, dim = 3) 
#' } 
#' @export 
correct_image_dim <- function(image, dim = 3){
  out.img <- image
  while (length(dim(out.img)) > dim) {
    out.img <- drop_img_dim(out.img)
  }
  return(out.img)
}