#' @title SPM12 Normalize (Estimate)
#'
#' @description Estimate SPM12 (Spatial) Normalization from image
#' @param filename File to be normalized to the template
#' @param add_spm_dir Add SPM12 directory from this package
#' @param spmdir SPM dir to add, will use package default directory 
#' @param clean Remove scripts from temporary directory after running
#' @param verbose Print diagnostic messages
#' @param ... Arguments passed to \code{\link{run_spm12_script}}
#' @export
#' @import fslr
#' @import matlabr
#' @return Result from run_matlab_script 
spm12_normalize_est <- function(filename, 
                                add_spm_dir = TRUE,
                                spmdir = spm_dir(),                          
                                clean = TRUE,
                                verbose = TRUE,
                                ...
){
  
  install_spm12()
  # check filenames
  filename = filename_check(filename)
  
  jobvec = c(filename, spmdir)
  names(jobvec) = c("%filename%", "%spmdir%")
  
  res = run_spm12_script( script_name = "Normalize_Estimate",
                          jobvec = jobvec,
                          mvec = NULL,
                          add_spm_dir = add_spm_dir,
                          spmdir = spmdir,
                          clean = clean,
                          verbose = verbose,
                          ...)
  return(res)
}


