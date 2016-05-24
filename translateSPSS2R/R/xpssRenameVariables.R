#' Renaming Variables
#'
#' R implementation of the SPSS \code{RENAME VARIABLES} function. xpssRenameVariables renames variables within an exisiting data.frame or xpssFrame object.
#'
#' Modifies names of one or more variables within a selected dataset. The arguments oldVarNames and newVarNames must have the same length.
#'
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @param oldVarNames atomic character or character vector with the names of the variables to rename. oldVarNames must be a variable that already exists in the data set. 
#' @param newVarNames atomic character or character vector with the new variable names. 
#' @return Returns the data with the renamed variables.
#' @author Andreas Wygrabek
#' @examples 
#' data(fromXPSS)
#' 
#' xpssRenameVariables(fromXPSS, 
#' oldVarNames= c("V1", "V2", "V3"), 
#' newVarNames= c("Manufacturer", " Car Type", "Country"))
#' @export
xpssRenameVariables <- function(x, oldVarNames = NULL, newVarNames = NULL){

  ####################################################################
  ####################################################################
  
  functiontype <- "DM"
  x <- applyMetaCheck(x)
  
  ####################################################################
  ####################################################################
  ####################################################################
  for(i in 1:length(oldVarNames)) {
    if(!(is.element(oldVarNames[[i]],names(x)))) {
      stop("The selected variable has to be in the dataset")
    }
  }
  
  oldVarNames <- c(oldVarNames)
  newVarNames <- c(newVarNames)
  colnames(x)[colnames(x) %in% oldVarNames] <- newVarNames
  
  pos <- which(colnames(x) %in% newVarNames)
  
  for(i in 1:length(newVarNames)) {
    attributes(x[[pos[i]]])$varname <- newVarNames[i]  
  }
  
  x <- applyAttributeDemerge(x)
  

  return(x)
}




