#' Save data to a NIfTI file
#' 
#' This function can be used to output the results of voxelwise RLRT or
#' smoothing.
#' 
#' 
#' @param arr a 3D or 4D array containing data to be saved.
#' @param name.nii filename, excluding the .nii extension.
#' @return None; a NIfTI file is created.
#' @author Lei Huang \email{huangracer@@gmail.com}, Philip Reiss
#' \email{phil.reiss@@nyumc.org} and Rong Jiao \email{jiaorong007@@gmail.com}
#' @seealso \code{\link{nii2R}}
#' @export
R2nii <-
function(arr, name.nii) {
	  ndim = length(dim(arr))
    arr[is.na(arr)] = 0
    if (ndim==3) {
    	  nim = array(0, attributes(arr)$dim.nii)
        nim[attributes(arr)$x.ind, attributes(arr)$y.ind, attributes(arr)$z.ind] = arr
    }
    else if (ndim==4) {
    	  nim = array(0, c(attributes(arr)$dim.nii, dim(arr)[4]))
        nim[attributes(arr)$x.ind, attributes(arr)$y.ind, attributes(arr)$z.ind, ] = arr
    }
    
    oro.nifti::writeNIfTI(nim, name.nii, gzipped=FALSE, verbose=TRUE)
}



