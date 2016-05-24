#' @title AC/PC Reorientation
#' @description Function that AC/PC re-orients the images for SPM 
#' spatial normalization routine.  Uses nii_setOrigin from 
#' http://www.mccauslandcenter.sc.edu/CRNL/sw/spm8/spm.zip
#' @param infiles (character) Files to reorient.  
#' First file will be used to 
#' estimate AC/PC, then rest will be transformed
#' @param modality T1, T2, CT, fMRI (T2*)
#' @param spmdir (character) path for SPM8.  If NULL, assumes 
#' SPM8 is in matlabpath and so is spm8/toolbox
#' Must have nii_setOrigin installed.  In 
#' \code{system.file("", package="cttools")} from
#' http://www.mccauslandcenter.sc.edu/CRNL/sw/spm8/spm.zip
#' @param verbose (logical) Print diagnostic output
#' @param ... Arguments to pass to \code{\link{run_matlab_code}}
#' @return Exit code from MATLAB.  If not zero, there was an error
#' @import fslr
#' @import stringr
#' @export
acpc_reorient <- function(
  infiles, 
  modality = c("T1", "T2", "CT", "fMRI"),
  spmdir = spm_dir(), 
  #   add_cttools = TRUE,
  verbose=TRUE,
  ...
){
  
  install_spm12()

  
  infiles = checknii(infiles)
  if (verbose) cat(paste0("\n #Reorientation ", infiles[1], "\n"))
  matcmd = get_matlab()
  ### gantry tilt correction - make new folder
  ### ranem old folder - zip it and then run matlab script
  
  modality = match.arg(modality, 
                       c("T1", "T2", "CT", "fMRI"))
  modality_num = switch(modality,
                        "T1"=1,
                        "T2"=2,
                        "CT"=3,
                        "fMRI"=4)
  
  
  cmd = NULL
  if (!is.null(spmdir)){
    spmdir = path.expand(spmdir)
    cmd <- paste(cmd, sprintf("addpath(genpath('%s'));", spmdir))
  }
  
  #   if (add_cttools){
  #     cmd <- paste(cmd, sprintf("addpath('%s');", 
  #                             system.file("", package="cttools")))  
  #   }
  #   cmd <- paste(cmd, sprintf("addpath('%s/toolbox/rorden');", spmdir))
  
  limgs = length(infiles)
  imgs = sprintf("'%s',", infiles[1])
  if (limgs > 1){
    for (ifile in 2:limgs){
      imgs = paste( imgs, sprintf("'%s',", 
                                  infiles[ifile]))
    }
  }
  
  imgs = str_trim(imgs)
  imgs = gsub(",$", "", imgs)
  cmd <- paste(cmd, sprintf("runimgs = strvcat(%s);", imgs))
  cmd <- paste(cmd, paste0("nii_setOrigin(runimgs, ", 
                           modality_num,
                           ");"))
  x = run_matlab_code(cmd, ...)
  return(x)
}
