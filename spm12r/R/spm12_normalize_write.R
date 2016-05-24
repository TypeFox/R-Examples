#' @title SPM12 Normalize (Write)
#'
#' @description Applies SPM12 (Spatial) Normalization to images
#' @param filename Filename of deformation (nifti)
#' @param other.files Files to be written using the estimated 
#' normalization
#' @param retimg Logical indicating if image should be returned or
#' result from \code{\link{run_matlab_script}}
#' @param reorient if \code{retimg=TRUE} pass to \code{\link{readNIfTI}}
#' @param add_spm_dir Add SPM12 directory from this package
#' @param spmdir SPM dir to add, will use package default directory 
#' @param clean Remove scripts from temporary directory after running
#' @param verbose Print diagnostic messages
#' @param ... Arguments passed to \code{\link{run_spm12_script}}
#' @export
#' @import fslr
#' @import matlabr
#' @return Result from run_matlab_script 
spm12_normalize_write <- function(filename, 
                            other.files = NULL,
                            retimg = TRUE,                            
                            reorient = FALSE,
                            add_spm_dir = TRUE,
                            spmdir = spm_dir(),                          
                            clean = TRUE,
                            verbose = TRUE,
                            ...
){
  install_spm12()
  
  # check filenames
  filename = filename_check(filename)
  
  if (!is.null(other.files)){
    other.files = filename_check(other.files)
  } else {
    stop("No files specified, no files written")
  }
  other.fnames = other.files
  # Pasting them together
  other.files = sapply(other.files, function(o){
    paste0("'", o, "'\n")
  })
  other.files = paste(other.files,  collapse = " ")
  
  jobvec = c(filename, other.files)
  names(jobvec) = c("%filename%", "%resample%")
  
  res = run_spm12_script( script_name = "Normalize_Write",
                          jobvec = jobvec,
                          mvec = NULL,
                          add_spm_dir = add_spm_dir,
                          spmdir = spmdir,
                          clean = clean,
                          verbose = verbose,                           
                          ...)
  
  if (retimg){
    outfiles = file.path(dirname(other.fnames), 
                         paste0("w", basename(other.fnames)))  
    res = lapply(outfiles, readNIfTI, reorient = reorient)
    return(res)
  }  
  return(res)
}


