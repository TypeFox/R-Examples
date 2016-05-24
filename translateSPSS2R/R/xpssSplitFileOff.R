#' Set the data split off
#'
#' R implementation of the SPSS \code{SPLIT FILE} argument
#'
#' @usage xpssSplitFileOff(x)
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @return The function return the input data with modified attributes.
#' @author Bastian Wiessner
#' @seealso \code{\link{xpssSplitFile}} \code{\link{xpssFilterOff}}
#' \code{\link{xpssTemporary}} 
#' @export
xpssSplitFileOff <- function(x){
  
  ###################################################################
  ####################################################################
  
  functiontype <- "DM"
  x <- applyMetaCheck(x)
  
  ####################################################################
  ####################################################################
  ####################################################################
  
  attributes(x)$SPLIT_FILE <- NULL
  attributes(x)$SPLIT_FILE <- F
  
  print("Data split is deactivated")
  return(x)
}


