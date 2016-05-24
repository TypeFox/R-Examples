#' @title SPM12 Smooth
#'
#' @description Performs SPM12 Smoothing on an Image
#' @param filename File to be smoothed
#' @param retimg Logical indicating if image should be returned or
#' result from \code{\link{run_matlab_script}}
#' @param fwhm Full-Width Half Max to smooth
#' @param prefix Prefix to append to front of image filename
#' @param add_spm_dir Add SPM12 directory from this package
#' @param spmdir SPM dir to add, will use package default directory 
#' @param clean Remove scripts from temporary directory after running
#' @param verbose Print diagnostic messages
#' @param reorient if \code{retimg=TRUE} pass to \code{\link{readNIfTI}}
#' @param ... Arguments passed to \code{\link{run_spm12_script}}
#' \code{\link{readNIfTI}}
#' @export
#' @import fslr
#' @import matlabr
#' @return Result from run_matlab_script or nifti file, depending on 
#' \code{retimg}
spm12_smooth <- function(filename, 
                         retimg = TRUE,
                         fwhm = 8,
                         prefix = "s",
                         add_spm_dir = TRUE,
                         spmdir = spm_dir(),                          
                         clean = TRUE,
                         verbose = TRUE,
                         reorient = FALSE,
                         ...
){
  install_spm12()  
  # check filenames
  filename = filename_check(filename)
  
  jobvec = c(filename, prefix, fwhm)
  names(jobvec) = c("%filename%", "%prefix%", "%fwhm%")
  
  res = run_spm12_script( script_name = "Smooth",
                          jobvec = jobvec,
                          mvec = NULL,
                          add_spm_dir = add_spm_dir,
                          spmdir = spmdir,
                          clean = clean,
                          verbose = verbose,
                          ...)
  outfile = file.path(dirname(filename), 
                      paste0(prefix, basename(filename)))  
  if (retimg){
    res = readNIfTI(outfile, reorient = reorient)
    return(res)
  }
  return(outfile)
}


