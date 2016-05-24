##' A function to standardize the unit
##'
##' The function standardize the data to the desirable unit when the
##' multiplier vector is supplied. For example per 1000 people is
##' scaled to per person by supplying a multiplier of 1000.
##'
##' @param df The data frame containing the data to be scale
##' @param multiplier The named vector with the multiplier to be
##' scaled. The name is mandatory in order for the function to identify
##' the variable in the data frame. A data.frame can also be supplied
##' with the first column being the name and the second being the
##' numeric multiplier.
##' @export
##' @examples
##'
##' ## Create the data frame
##' test.df = data.frame(FAOST_CODE = 1:5, Year = 1995:1999,
##'   var1 = 1:5, var2 = 5:1)
##'
##' ## Create the named vector for scaling
##' multiplier = c(1, 10)
##' names(multiplier) = c("var1", "var2")
##'
##' ## Scale the data
##' scaleUnit(test.df, multiplier = multiplier)

scaleUnit = function(df, multiplier){
  printLab("Scaling to base unit")
  cat("\t\tNOTE: This function should only be performed once\n\n")
  if(is.data.frame(multiplier) & NCOL(multiplier) == 2){
    colModes = sapply(multiplier, mode)
    if(any(colModes == "numeric") && any(colModes == "character")){
      multiplierVec = multiplier[, which(colModes == "numeric")]
    names(multiplierVec) = multiplier[, which(colModes == "character")]
    } else {
      stop("Invalid multiplier data.frame")
    }
  } else if(is.vector(multiplier)){
    multiplierVec = multiplier
  } else {
    stop("Invalid multiplier input")
  }
  if(is.null(names(multiplierVec)))
    stop("Please supply names to the multiplierVec vector")
  multiplierVec[is.na(multiplierVec)] = 1
  inIndex = which(names(multiplierVec) %in% colnames(df))
  if(length(inIndex) != length(multiplierVec)){
    cat("The following variables are not in the data frame:\n")
    cat(paste("\t", names(multiplierVec[-inIndex]), "\n", sep = ""))
  }
  n = NROW(df)
  multiplierVec.mat = matrix(rep(multiplierVec[inIndex], n), nrow = n,
    byrow = TRUE)
  df[, names(multiplierVec[inIndex])] = df[, names(multiplierVec[inIndex])] *
    multiplierVec.mat
  df
}
