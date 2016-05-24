#' Timer function (as in MATLAB)
#' 
#' Funtion to stop a timer.  Start with tic().
#' 
#' @param echo Print time to screen?
#' @param name The saved name of the time object.
#' 
#' @note This is a modified version of the same function in the matlab package \code{\link[matlab]{toc}}
#' 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_toc.R
#' @export
## Function written to match MATLAB function
## Author: Andrew Hooker

#toc <- function(...){
#    return(matlab::toc(echo=FALSE))
#}

toc <- function (echo = TRUE,name=".poped_savedTime") 
{
  #prevTime <- get(name, pos = 1)
  prevTime <- get(name, envir=.PopedNamespaceEnv)
  
  diffTimeSecs <- proc.time()[3] - prevTime
  if (echo) {
    cat("Elapsed time:", diffTimeSecs, "seconds.\n") 
    #cat(sprintf("Elapsed time is %f seconds.\n", diffTimeSecs)) 
    #cat(sprintf("Elapsed time is %f seconds.", diffTimeSecs),"\n")
    return(invisible(diffTimeSecs))
  }
  else {
    return(diffTimeSecs)
  }
}
