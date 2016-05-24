#' Modifies value labels
#'
#' R implementation of the SPSS \code{VALUE LABEL} function. xpssValueLabels creates value labels for specific variables. The values of the label get stored in attributes of the variable.
#'
#' @param x a (non-empty) data.frame or input data of class "xpssFrame". 
#' @param variables atomic character or character vector with the names of the variables.
#' @param values atomic numeric or numeric vector containing the values of the variable.
#' @param labels atomic character or character vector containing the variable labels.
#' @param datevariables atomic date or date vector with the names of the date variables.
#' @param datevalues atomic date or date vector containing the value of the date, the values has to be like the old date format.
#' @param datelabels atomic character or character vector containing the date labels.
#' @details The SPSS variables are stored at the variable itself. 
#' \cr In contrast to \code{\link{xpssAddValueLabels}} , \code{\link{xpssValueLabels}} does erase existing value labels. \cr If the value label for a specific variable already exists, all value labels for that variable get overwritten. \cr If the value label for a specific variable does not exist, the value label gets created and all existing value labels for that variable get deleted. 
#' \cr\cr A variable can have the following attributes: 
#' \code{value.labels}, \code{defined.MIS}, \code{MIS}, \code{varname}, \code{variable.label}
#' @author Andreas Wygrabek
#' @seealso \code{\link{xpssAddValueLabels}} \code{\link{xpssVariableLabels}}
#' @examples
#' 
#' data(fromXPSS)
#' 
#' temp <- xpssValueLabels(fromXPSS, 
#'                            variables = "V1", 
#'                            value = 1 ,
#'                            label = "Label1")
#'
#' fromXPSS <- xpssAddValueLabels(fromXPSS, 
#'                            variables = "V1", 
#'                            values = c("A","B"),
#'                            labels = c("Label1","Label2"))                           
#'                            
#' @export
xpssValueLabels <- function(x, variables = NULL, values = NULL, labels = NULL,
                            datevariables = NULL, datevalues = NULL, datelabels = NULL){
    
  ####################################################################
  ####################################################################
  
  functiontype <- "DM"
  x <- applyMetaCheck(x)
  
  ####################################################################
  ####################################################################
  ####################################################################
  for(i in 1:length(variables)) {
    if(length(variables>0) && (!(is.element(variables[[i]],names(x))))) {
      stop("The selected variables has to be in the dataset")
    }  
  }
  for(i in 1:length(datevariables)) {
    if(length(datevariables>0) && (!(is.element(datevariables[[i]],names(x))))) {
      stop("The selected datevariables has to be in the dataset")
    }  
  }
  
    
    if(length(values) != length(labels)){
        stop("Vectors values and labels dont have the same length")
    }
    
    if("NA" %in% labels){
        stop("Using NA as label is not possible")
    }
    if(!is.null(values))
    {
      if(class(values) != "numeric" && class(values) != "character"){
        stop("Only numeric or character values are allowed for the argument values")
      }
    } 
   if(length(variables) > 0) {
       for(i in 1:length(variables)){
         
         if(class(variables[i]) != "numeric" | class(variables[i]) != "factor"){
           "Input-Variables from variables have to be numeric or factor"
         }
         
         names(values) <- labels
         attr(x[,variables[i]], "value.labels") <- values
       }
     }
     
     if(length(datevariables) >0 ){
       for(i in 1:length(datevariables)){
         
         if(class(datevariables[i]) != "date" | class(datevariables[i]) != "POSIXlt" | class(datevariables[i]) != "POSIXt" | class(datevariables[i]) != "POSIXct"){
           "Input-Variables from datevariables have to be class date or POSIX"
         }
         
         names(datevalues) <- datelabels
         attr(x[,datevariables[i]], "value.labels") <- datevalues
       }
     }
   
   
   x <- applyAttributeDemerge(x)
   
   return(x)
   }
    



    
    
    
