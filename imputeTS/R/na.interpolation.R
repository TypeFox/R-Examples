#' @title Missing Value Imputation by Interpolation
#' 
#' @description Uses either linear or spline interpolation to replace missing values.
#'  
#' @param x Numeric Vector (\code{\link{vector}}) or Time Series (\code{\link{ts}}) object in which missing values shall be replaced
#' @param option Algorithm to be used. Accepts the following input:
#' \itemize{
#'    \item{"linear" - for linear interpolation using \link{approx} }
#'    \item{"spline" - for spline interpolation using \link{spline}}
#'    \item{"stine" - for Stineman interpolation using \link{stinterp}}
#'    }
#' @param ... Additional parameters to be passed through to \link{approx} or \link{spline} interpolation functions
#' @return Vector (\code{\link{vector}}) or Time Series (\code{\link{ts}}) object (dependent on given input at parameter x)
#' 
#' @details Missing values get replaced by values of a \link{approx} or \link{spline} interpolation.
#'  Both functions are thereby called with their default parameters (except for parameter rule at approx, which is set to 2)
#'  If more flexibility in parameter settings is needed, na.approx / na.spline from the \pkg{zoo} 
#'  package can be a option.
#' 
#' @author Steffen Moritz
#' @seealso  \code{\link[imputeTS]{na.kalman}}, \code{\link[imputeTS]{na.locf}},
#'  \code{\link[imputeTS]{na.ma}}, \code{\link[imputeTS]{na.mean}},
#'  \code{\link[imputeTS]{na.random}}, \code{\link[imputeTS]{na.replace}},
#'  \code{\link[imputeTS]{na.seadec}}, \code{\link[imputeTS]{na.seasplit}}
#'  
#' @examples
#' #Prerequisite: Create Time series with missing values
#' x <- ts(c(2,3,4,5,6,NA,7,8))
#' 
#' #Example 1: Perform linear interpolation
#' na.interpolation(x)
#' 
#' #Example 2: Perform spline interpolation
#' na.interpolation(x, option ="spline")
#' 
#' #Example 3: Perform stine interpolation
#' na.interpolation(x, option ="stine")
#' 
#' @import stats
#' @import stinepack
#' @export


na.interpolation <- function(x, option = "linear", ...) { 
 
  data <- x
  
  #Check for wrong input 
  data <- precheck(data)
  
  #if no missing data, do nothing
  if(!anyNA(data)) {
    return(data)
  }
  
  ##
  ## Imputation Code
  ##
  missindx <- is.na(data)  

  n <- length(data)
  
  allindx <- 1:n
  indx <- allindx[!missindx]
  
  if(option =="linear") {
    interp <- as.ts(approx(indx,data[indx],1:n, rule=2, ...)$y)
  }
  else if(option == "spline") {
    interp <- as.ts(spline(indx,data[indx],n = n, ... )$y)
  }
  else if(option == "stine") {
    interp <- as.ts(stinterp(indx,data[indx],1:n, ...)$y)
    #avoid NAs at the beginning and end of series // same behavior like for approx with rule = 2.
    if(any(is.na(interp))) {interp <- na.locf(interp, na.remaining= "rev")}
  }
  else {
    stop("Wrong parameter 'option' given. Value must be either 'linear' or 'spline'.")
  }
  
  for (i in 1:length(data)) {
    if (is.na(data[i])) {
      data[i] <- interp[i]
    }
  }
  
  return(data)
}