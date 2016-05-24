#' @title SPM12 Coregister (Estimate and Reslice)
#'
#' @description Performs SPM12 coregistration estimation and reslicing on an Image
#' @param fixed File that is assumed fixed
#' @param moving moving file to be registered to fixed space
#' @param other.files Other files to register to fixed, in same space as moving file
#' @param prefix Prefix to append to front of image filename 
#' @param add_spm_dir Add SPM12 directory from this package
#' @param spmdir SPM dir to add, will use package default directory 
#' @param clean Remove scripts from temporary directory after running
#' @param verbose Print diagnostic messages
#' @param outdir Directory to copy results.  If full filename given, then results will
#' be in \code{dirname(filename)}
#' @param ... Arguments passed to \code{\link{run_spm12_script}}
#' @export
#' @import fslr
#' @import matlabr
#' @return Result from run_matlab_script 
spm12_coregister <- function(fixed,
                             moving, 
                             other.files = NULL,
                             prefix = "r",
                             add_spm_dir = TRUE,
                             spmdir = spm_dir(),                          
                             clean = TRUE,
                             verbose = TRUE,
                             outdir = NULL,                       
                             ...
){
  
  install_spm12()
  
  ########################
  # Getting Number of Time points
  ########################  
  
  # check filenames
  fixed = filename_check(fixed)
  moving = filename_check(moving)
  
  if (is.null(other.files)){
    other.files = "''"
    other.ofiles = NULL
    other = FALSE
  } else {
    other.files = filename_check(other.files)
    other.ofiles = file.path(dirname(other.files),
                             paste0(prefix, basename(other.files)))    
    other.files = rvec_to_matlabcell(other.files)
    other = TRUE
  }
  omoving = file.path(dirname(moving),
                      paste0(prefix, basename(moving)))      
  ##########################################################
  # Pasting together for a 4D file
  ##########################################################
  
  
  jobvec = c(fixed, moving, other.files, prefix)
  names(jobvec) = c("%reffile%", "%sourcefile%", "%otherfile%", "%prefix%")  
  
  res = run_spm12_script( script_name = "Coregister",
                          jobvec = jobvec,
                          mvec = NULL,
                          add_spm_dir = add_spm_dir,
                          spmdir = spmdir,
                          clean = clean, 
                          verbose = verbose, 
                          ...)
  stopifnot(res == 0)
  ####################
  # Copy outfiles
  ####################  
  if (!is.null(outdir)){
    file.copy(omoving, to = outdir, overwrite = TRUE)
    if (other) {
      file.copy(other.ofiles, to = outdir, overwrite = TRUE)
    }
  }
  
  l = list(outfile = moving, 
           other.outfiles = other.ofiles)
  return(l)
}


