#' Standardize Variables
#'
#' This function takes in a dataframe, the name of any number of variables.  It then returns either a vector or a dataframe of scaling results.  If passed a single variable name, standardize will return a the standardized variable as a vector, when passed 2 or more variable names, standardize will return a data frame containing all of the standardized variables. 
#'
#' Often times we are forced to compare variables which exist on different scales, but how do you compare the coefficient for a country's population to one that's much smaller?  Standardizing your variables can make reading regression results more useful because it will make coeficients more comparable.
#'
#' @param data data.frame object that contains both the dependent variable and predictor variables you want to regress.
#' 
#' @param variable variable name or vector of names for variables you want standardized.
#' 
#' @param type by default, 'absolute' will scale every variable from 0 to 1.  "classic" is a little more complicated where the variable is rescaled the mean equaling 0 and a standard deviation is 1.
#' 
#' @return if the function is passed a single variable to standardize, it will return a vector of all obeservations in the variable standardized.  If the function is passed more than one variable name, it will return a data-frame containing the new observation values.  
#'
#' @export
#' 
#' @examples 
#' exam.df <- iris
#' standardize(exam.df, "Petal.Width")
#' standardize(exam.df, c("Petal.Width", "Petal.Length"), type = "classic")

standardize <- function(data, variable, type = "absolute") {
  
  if(!is.character(variable)) {
    warning("variable must be entered as character object")
  }
  if(length(variable) > 1) {
    if(type == "absolute") {
      createdvars <- list()
      for(i in variable){
        name <- paste("stand.", i, sep = "")
        createdvars[[length(createdvars)+1]] <- name
        assign(name, ( data[,i] - min(data[,i]) ) / ( max(data[,i]) - min(data[,i]) ))
      }
    }
    else{if(type == "classic") {
      createdvars <- list()
      for(i in variable){
        name <- paste("stand.", i, sep = "")
        createdvars[[length(createdvars)+1]] <- name
        assign(name, ( data[,i] - mean(data[,i]) ) / sd(data[,i]))
        }
      }
    }
    stand.vars<- data.frame(data[,variable])
    for(i in createdvars){
      stand.vars <- cbind(stand.vars, get(i))
    }
    for(i in variable){
      stand.vars[,i] <- NULL
      }
    colnames(stand.vars) <- createdvars
    return(stand.vars)
  }
  if(length(variable) == 1)
    if(type == "absolute") {
      temp <- ( data[,variable] - min(data[,variable]) ) / ( max(data[,variable]) - min(data[,variable]) )
    }
    
    if(type == "classic") {
      temp <- ( data[,variable] - mean(data[,variable]) ) / sd(data[,variable]) 
    }
  return(temp)
}
