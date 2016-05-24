#' Timer function (as in MATLAB)
#' 
#' Funtion to start a timer.  Stop with toc().
#' 
#' @param gcFirst Perform garbage collection?
#' @param name The saved name of the time object.
#' 
#' @note This is a modified version of the same function in the matlab package \code{\link[matlab]{tic}}
#' 
#' @family MATLAB
#' @example tests/testthat/examples_fcn_doc/examples_tic.R
#' @export
## Function written to match MATLAB function
## Author: Andrew Hooker

#tic <- function(...){
#    return(matlab::tic())
#}

tic <- function (gcFirst = FALSE,name=".poped_savedTime") 
{
  if (gcFirst == TRUE) {
    gc(verbose = FALSE)
  }
  
  #   if(exists(envir)){
  #     if(!is.environment(eval(parse(text=envir)))) eval(parse(text=envir)) <- new.env(parent=baseenv())
  #   } else {
  #.PopEDNamespaceEnv <- new.env(parent=.GlobalEnv)
  #   }
  #   if(!exists(name,envir=.PopEDNamespaceEnv)){
  #   .PopEDNamespaceEnv <- new.env(parent=baseenv())
  #   }
  assign(name, proc.time()[3], envir = .PopedNamespaceEnv)
  
  #assign(name, proc.time()[3], pos = package:PopED)
  
  #assign(name, proc.time()[3], pos = eval(parse(text="1")))
  invisible()
}