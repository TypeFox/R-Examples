#' @title Perform filename checks for SPM
#'
#' @description Checks a filename to see if nifti and expands paths
#' to absolute paths
#' @param filename filename of an image or nifti object
#' @param ... arguments passed to \code{\link{checknii}}
#' @export
#' @import fslr
#' @seealso \code{\link{checknii}}
#' @return Character of filename
filename_check <- function(filename, # filename of an image
                           ... # arguments passed to \code{\link{checknii}}
                           ){
  ###########################
  # Passing to see if image or filename passed in
  ###########################  
  filename = checknii(filename, ...)
  filename = path.expand(filename)
  filename = normalizePath(filename)

  stopifnot(inherits(filename, "character"))
  stopifnot(all(file.exists(filename)))
  return(filename)
}