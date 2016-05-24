#' Creates a FILTER subset 
#'
#' R implementation of the SPSS \code{FILTER} argument.
#'
#' \code{xpssFilter} creates a subset of the actual dataset, without deleting the excluded variables, respectively without deleting the excluded values. After activating \code{xpssFilter}, only the subset will be used, the excluded data will be ignored for the following actions until  \code{\link{xpssFilterOff}} terminates the filtering. As noted before those cases are \strong{not actually} deleted and will be available after the filter is turned off. 
#' \cr\cr \strong{Important:} 
#' \cr All changes are used on the complete dataset, except for the function being an \emph{data exploring} or \emph{data analyzing} function \cr \cr
#'  \tabular{rlll}{
#'  \tab Type of Function \tab Example Function \tab Dataset Usage \cr
#' \tab Data Management  \tab \code{\link{xpssSelectIf}} \tab Uses thecomplete dataset\cr
#' \tab Data Modifing \tab \code{\link{xpssRecode}} \tab Uses the complete dataset\cr
#' \tab Data Exploring \tab \code{\link{xpssDescriptives}} \tab Uses theworking dataset only\cr
#' \tab Data Analyzing \tab \code{\link{xpssRegression}} \tab Uses the working dataset only\cr
#'}
#' @usage xpssFilter(x,variable = NULL, filtervalue = 1)
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @param variable atomic character with the name of the variable. 
#' @param filtervalue atomic character or atomic numeric which contains the filtervalue.
#' @return Output is a subset of the actual dataset under the predetermined condition of the filtervalue.
#' @author Bastian Wiessner
#' @seealso Related Functions \code{\link{xpssDoIf}} \code{\link{xpssFilterOff}} \code{\link{xpssSample}} \code{\link{xpssSelectIf}} \code{\link{xpssTemporary}} 
#' @examples 
#' data(fromXPSS)
#' 
#' fromXPSS <- xpssFilter(x=fromXPSS, variable = "V3", filtervalue=1)
#' 
#' xpssDescriptives(x=fromXPSS, variables = "V6")
#' 
#' xpssFilterOff(x=fromXPSS)
#' @export
xpssFilter <- function(x, variable = NULL, filtervalue = 1){

  if(!(is.element(variable,names(x)))) {
    stop("The selected variable has to be in the dataset")
  }
  
  functiontype <- "ME"
  x <- applyMetaCheck(x)
    
  #Exception if the old Filter gets replaced by a new Filter
  if(is.null(attributes(x)$FILTERED_DATA) == FALSE) {
    x <- xpssFilterOff(x)   
  }  
  pos <-  which(x[,variable] == filtervalue)
  
  attr(x, "FILTER") <- paste0(variable, " == ", filtervalue)
  attr(x, "FILTERED_DATA") <- x[setdiff(1:nrow(x),pos),]
  
  attBack <- attributesBackup(x)
  
  x <- x[pos,]
  
  print("Filter is activated: Keep the rownames to switch the filter off")
  
  x <- applyAttributes(x, attBack)
  
  return(x)
}




