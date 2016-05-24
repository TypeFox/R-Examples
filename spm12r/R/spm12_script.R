#' @title Find SPM12 Script
#'
#' @description Copies the SPM12 script from the scripts directory
#' to a temporary file
#' @param script_name Name of the script filename without ".m" ext
#' @param outdir Path to copy scripts and run
#' @export
#' @return Chracter vector of script paths
#' @examples spm12_script(script_name = "Segment")
spm12_script <- function(script_name, outdir = tempdir()){

    m_scripts = system.file("scripts", 
                           paste0(script_name, c(".m")), 
                           package="spm12r")
  
  ####################
  # Use General Executable
  ####################  
  miss = m_scripts %in% ""
  if (any(miss)){
    m_scripts[miss] = system.file("scripts", "Executable.m", 
                    package="spm12r")
  }
  ####################
  # Get Jobfile
  ####################    
  job_scripts = system.file("scripts", 
                        paste0(script_name, c("_job.m")), 
                        package="spm12r")  
  scripts = c(job = job_scripts, script = m_scripts)
  scripts = scripts[scripts != "", drop = FALSE]
  nn = names(scripts)
  if (length(scripts) > 0){
    file.copy(scripts, to = outdir, overwrite = TRUE)
    scripts = file.path(outdir, basename(scripts))
    names(scripts) = nn
  }
  return(scripts)
}
