#' Warning function
#' 
#' Writes or stores messages for various situations that can occur
#' 
#' @param warning.msg the message to be stored/printed (optional)
#' @param warning.mode report or print errors (default report)
#' @param MAE.warnings character vector of existing warning messages
#' @return None
#' @author Dave Miller & Laura Marshall
mae.warning <- function(warning.msg=NULL, warning.mode="store", MAE.warnings)
#
# errors
#
# A nice, neat way of storing and printing errors/warnings.
#
# Arguments:
# 
#  errmode	- report or print errors (default report)
#  errmsg	- the message to be stored/printed (optional)
#
# Values:
#
#  We return either a formatted error message, TRUE or a list()
#  of all the errors.
#
# dlm 25-Aug-05  Initial work started. At the moment we are just
#		 able to format the data, nothing more complicated
#		 yet.
#		 By default in future it should store the errors in
#		 a global variable.
# lhm 4-Apr-12 Now storing in a global variable
#
{

  if(warning.mode == "report"){
    warning(warning.msg)
  }else if(warning.mode == "store"){
    MAE.warnings <- c(MAE.warnings, warning.msg)  
  }else{
    warning("\nOnly report and store modes are implemented.\n")
  }
  return(MAE.warnings)     
}
