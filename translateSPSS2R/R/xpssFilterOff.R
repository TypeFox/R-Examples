#' Ends a FILTER subset 
#'
#' R implementation of the SPSS \code{FILTER OFF} Function.
#'
#'  xpssFilterOff terminates the filtering and merges the excluded data with the actual subset of the dataset.
#' \cr\cr \strong{Important:} 
#' \cr All changes are used on the complete dataset, except for the function beeing an \emph{data exploring} or \emph{data analyzing} function. \cr \cr
#'  \tabular{rlll}{
#'  \tab Type of Function \tab Example Function \tab Dataset Usage \cr
#' \tab Data Management  \tab \code{\link{xpssSelectIf}} \tab Uses the complete dataset\cr
#' \tab Data Modifing \tab \code{\link{xpssRecode}} \tab Uses the complete dataset\cr
#' \tab Data Exploring \tab \code{\link{xpssDescriptives}} \tab Uses the working dataset only\cr
#' \tab Data Analyzing \tab \code{\link{xpssRegression}} \tab Uses the working dataset only\cr
#'}
#' \strong{NOTE:} For temporary case selection, specify \code{xpssTemporary} before \code{xpssDoIf}.
#'
#' @usage xpssFilterOff(x)
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @return Output is the original dataset.
#' @author Andreas Wygrabek
#' @examples
#' data(fromXPSS)
#' 
#' temp <- xpssDoIf(x=fromXPSS, cond = "V3 == 1")
#' 
#' temp <- xpssRecode(x=temp,variables="V5",rec="lo:78 = 1; else = 2")
#' 
#' temp <- xpssEndIf(x=temp)
#' 
#' @export
xpssFilterOff <- function(x){
  
  functiontype <- "ME"
  x <- applyMetaCheck(x)
    
  if(attributes(x)$FILTER == FALSE){
      stop("Filter not activated")
  }
  
  x <- rbind(x,attributes(x)$FILTERED_DATA)
  attBack <- attributesBackup(x)

  
  x <- x[order(as.numeric(rownames(x))),]  
  x <- applyAttributes(x, attBack)
  
  attributes(x)$FILTER  <-  FALSE
  attributes(x)$FILTERED_DATA <- NULL
  
  return(x)
}











