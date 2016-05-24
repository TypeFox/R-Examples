#' Summarises warnings
#'
#' Sumarises warnings generated during the bootstrap and removes the 
#' MAE.warnings global object.  
#' 
#' @param MAE.warnings character vector of warning messages 
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#'
process.warnings <- function(MAE.warnings){
# process.warnings function to summarise warnings
#
# Arguments: none
# Value: none
# Function Calls: none
#
  mae.warning <- unique(MAE.warnings)
  no.warnings <- NULL
  if(length(mae.warning) > 0){
    cat("Warning messages: \n")  
  }
  for(i in seq(along = mae.warning)){
    no.warnings <- length(which(MAE.warnings == mae.warning[i])) 
    cat(paste(i, ": ", mae.warning[i], " [warning occured ",no.warnings," times]\n", sep = "")) 
  }
}  


