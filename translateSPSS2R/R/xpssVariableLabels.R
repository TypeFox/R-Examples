#' Modifies variable labels
#'
#' R implementation of the SPSS \code{VARIABLE LABEL} function. Changing the label of a variable. In the structure of xpss-data the variable label is an attribute of each variable.
#'
#' @param x a (non-empty) data.frame or input data of class "xpssFrame". 
#' @param variables atomic character or character vector with the names of the variable(s).
#' @param labels atomic character of character vector with labels for the specified variables in variables. The labels are associated in order of appearence
#' of the variables.   
#' @return Input Data with modified attribute variable label.
#' @author Andreas Wygrabek
#' @seealso \code{\link{attributes}} \code{\link{attr}} \code{\link{xpssValueLabels}}
#' @examples
#' data(fromXPSS)
#' 
#' daten <- xpssVariableLabels(fromXPSS, c("V4", "V7_1"), c("Label1", "Label2"))
#' @export
xpssVariableLabels <- function(x,  variables = NULL, labels = NULL){
    
  ####################################################################
  ####################################################################
  
  functiontype <- "DM"
  x <- applyMetaCheck(x)
  
  ####################################################################
  ####################################################################
  ####################################################################
  
  
    myList <- as.list( variables)
  
    if(any(!(is.element(variables,names(x))))) {
      stop("The selected variable has to be in the dataset")
    }
    
    if(length( variables) != length(labels)){
        stop("Length of  variables and labels has to be the same")
    }
    
    for(i in 1:length( variables)){
    attr(x[, variables[[i]]], "variable.label") <- labels[[i]]
    }
  
  
  x <- applyAttributeDemerge(x)
  
    return(x)
}
    




