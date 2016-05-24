#' @title SPM12 Realign (Estimate)
#'
#' @description Performs SPM12 Realignment estimation on an Image
#' @param filename Files to be realigned 
#' @param fwhm Full-Width Half Max to smooth 
#' @param register_to Should the files be registered to the first or the mean
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
#' @return Character list of filenames from output 
spm12_realign_est <- function(filename, 
                              fwhm = 5,                              
                              register_to = c("first", "mean"),
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
  if (verbose){
    cat("# Getting Number of Time Points\n")
  }
  time_points = ntime_points(filename)
  
  # check filenames
  filename = filename_check(filename)
  stub = nii.stub(filename, bn=TRUE)[1]
  rpfile = file.path(dirname(filename),
                     paste0("rp_", stub, ".txt"))
  #   meanfile = file.path(dirname(filename),
  #                        paste0("mean", stub, ".nii"))
  matfile = file.path(dirname(filename),
                      paste0(stub, ".mat"))
  
  ##########################################################
  # Pasting together for a 4D file
  ##########################################################
  filename = paste0(filename, ",", time_points)
  filename = rvec_to_matlabcell(filename)
  
  register_to = match.arg(register_to, c("first", "mean"))
  register_to = switch(register_to,
                       first = 0, 
                       mean = 1)
  
  jobvec = c(filename, fwhm, 
             register_to, spmdir)
  names(jobvec) = c("%filename%", "%fwhm%", 
                    "%registerto%", "%spmdir%")  
  
  res = run_spm12_script( script_name = "Realign_Estimate",
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
    file.copy(rpfile, to = outdir, overwrite = TRUE)
    #     file.copy(meanfile, to = outdir, overwrite = TRUE)
    file.copy(matfile, to = outdir, overwrite = TRUE)    
  }
  outfiles = c(rp = rpfile, 
               #                mean = meanfile, 
               mat = matfile)
  return(outfiles)
}


