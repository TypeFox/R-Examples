##' This is a function that use to match original data and addin data with identified variable.
##'
##' This function is for internal use only, so no need to export it. It figures out the ID in the "addin" file then merge variables in addin file to the original file. This function is used for further "lapply" porcess.
##'
##' @title An internal function for data merging.
##'
##' @param row One sample that is already divided from the original file.
##' @param data The "addin" file.
##' @param mergeVar The variable that use to merge.
##'
##' @return is single line contains original variables and addin variables.
##' @author Xiaorui.Zhu
consolida <- function(row, data, mergeVar){
  eachID <- row[mergeVar]
  pos <- match(eachID[1,], data[mergeVar][,1])
  return(cbind(row, data[pos,]))
}
