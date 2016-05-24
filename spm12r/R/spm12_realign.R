#' @title SPM12 Realign (Estimate and Reslice)
#'
#' @description Performs SPM12 realignment estimation and reslicing on an Image
#' @param filename Files to be realigned and resliced
#' @param fwhm Full-Width Half Max to smooth 
#' @param register_to Should the files be registered to the first or the mean
#' @param reslice Options for reslicing all - all images in filename,
#' 2:n - all images in filename 2:length(filename),
#' all+mean - all images and the mean, mean - mean only
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
spm12_realign <- function(filename, 
                          fwhm = 5,                              
                          register_to = c("first", "mean"),
                          reslice = c("all","2:n", "all+mean", "mean"),
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
  if (verbose){
    cat("# Getting Number of Time Points\n")
  }
  time_points = ntime_points(filename)
  
  reslice = match.arg(reslice, c("all","2:n", "all+mean", "mean"))
  
  # check filenames
  filename = filename_check(filename)
  ###################
  # If reslice is just mean, then the file is simply returned
  ###################  
  if (verbose){
    cat(paste0("# Reslice is ", reslice, "\n"))
  }
  if ( (reslice %in% "mean") ){
    outfile = filename
  } else {
    outfile = file.path(dirname(filename),
                        paste0(prefix, basename(filename)))    
  }
  
  reslice = switch(reslice,
                   "all" = "[2 0]", 
                   "2:n" = "[1 0]", 
                   "all+mean" = "[2 1]", 
                   "mean" = "[0 1]")
  
  stub = nii.stub(filename, bn=TRUE)[1]
  rpfile = file.path(dirname(filename),
                     paste0("rp_", stub, ".txt"))
  meanfile = file.path(dirname(filename),
                       paste0("mean", stub, ".nii"))
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

  jobvec = c(filename, prefix, fwhm, reslice,
             register_to, spmdir)
  names(jobvec) = c("%filename%", "%prefix%", "%fwhm%", "%reslice%",
                    "%registerto%", "%spmdir%")  
  
  res = run_spm12_script( script_name = "Realign",
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
    file.copy(outfile, to = outdir, overwrite = TRUE)
    file.copy(rpfile, to = outdir, overwrite = TRUE)
    if (!is.null(meanfile)){
      file.copy(meanfile, to = outdir, overwrite = TRUE)
    }
    file.copy(matfile, to = outdir, overwrite = TRUE)    
  }
  
  l = list(outfiles = outfile, 
           rp = rpfile, 
           mean = meanfile, 
           mat = matfile)
  return(l)
}


