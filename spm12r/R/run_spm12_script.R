#' @title Wrapper for running \code{spm12_script}
#'
#' @description Runs \code{\link{spm12_script}} with wrapper for 
#' spm12r functions
#' @param script_name Name of the script filename without .m ext,
#' passed to \code{\link{spm12_script}}
#' @param jobvec Vector of characters to be substituted in _job.m file
#' @param mvec Vector of characters to be substituted in .m file
#' @param add_spm_dir Add SPM12 directory from this package
#' @param spmdir SPM dir to add, will use package default directory 
#' @param clean Remove scripts from temporary directory after running
#' @param verbose Print diagnostic messages
#' @param ... Arguments to pass to \code{\link{spm12_script}}
#' @export
#' @return Result of \code{\link{run_matlab_script}}
run_spm12_script <- function(script_name, 
                            jobvec = NULL, 
                            mvec = NULL,
                            add_spm_dir = TRUE,
                            spmdir = spm_dir(),
                            clean = TRUE,
                            verbose = TRUE,
                            ...
                            ){
  install_spm12()
  
  scripts = spm12_script(script_name, ...)
  # put in the correct filenames
  job = readLines(scripts['job'])
  njvec = names(jobvec)
  if (any(is.na(jobvec))){
    print(jobvec)
    stop("There is an NA in jobvec")
  }
  for (ijob in seq_along(jobvec)){
    job = gsub(njvec[ijob], jobvec[ijob], job)
  }
  
  m = readLines(scripts['script'])  
  nmvec = names(mvec)
  for (ijob in seq_along(mvec)){
    m = gsub(nmvec[ijob], mvec[ijob], job)
  }
  m = gsub("%jobfile%", scripts['job'], m)
  
  #####################################
  # Adding in SPMDIR
  #####################################  
  if (add_spm_dir){
    if (verbose){
      cat(paste0("# Adding SPMDIR: ", spmdir, "\n"))
    }
    m = c(paste0("addpath(genpath('", spmdir, "'));"),
          m)
  }
  
  #####################################
  # Write Out files
  ##################################### 
  writeLines(m, con=scripts['script'])
  writeLines(job, con=scripts['job'])
  if (verbose){
    cat(paste0("# Running script", scripts['script'], "\nwhich calls ",
               scripts['job'], "\n"))
  }
  res = run_matlab_script(scripts['script'])
  if (verbose){
    cat(paste0("# Result is ", res, "\n"))
  }  
  #####################################
  # Cleaning up files
  ##################################### 
  if (clean) {
    file.remove(scripts)
    if (verbose){
      cat(paste0("# Removing scripts\n"))
    }      
  }
  return(res)
}