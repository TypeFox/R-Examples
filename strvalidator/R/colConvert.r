################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 06.10.2015: First version.

#' @title Convert Columns
#'
#' @description
#' Internal helper function.
#'
#' @details
#' Takes a data frame as input and return it after converting known numeric
#' columns to numeric.
#' 
#' @param data data.frame.
#' @param columns character string containing a regular expression
#'  (or character string for fixed = TRUE) to be matched in the given
#'  character vector (separate multiple column names by | in reg.exp).
#' @param ignore.case logical TRUE to ignore case in matching.
#' @param fixed logical TRUE if columns is a string to be matched as is.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame.
#' 
#' @export
#' 

colConvert <- function(data, columns="Height|Size|Data.Point",
                       ignore.case=TRUE, fixed=FALSE, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("data")
    print(str(data))
    print("columns")
    print(columns)
  }
  
  # Get all known numeric columns.
  selected <- grep(columns, names(data), ignore.case=ignore.case, fixed=fixed) 

  # Loop over all columns to change to numeric.
  for(c in seq(along=selected)){
    
    # Convert to numeric (handles factors using as.character).
    data[ , selected[c]] <- as.numeric(as.character(data[ , selected[c]]))
    
  }
  
  # Add attributes.
  attr(data, which="colConvert, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(data, which="colConvert, call") <- match.call()

  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }

  return(data)

}