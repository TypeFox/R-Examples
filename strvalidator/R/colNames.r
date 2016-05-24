################################################################################
# TODO LIST
# TODO: Make base name defining string optional instead of always "."

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom.
# 31.05.2015: New option 'numbered' for more precise control.
# 23.05.2015: Changed parameter name 'df' to 'data' for consistency with other functions.
# 15.12.2013: Fixed cropping column names with multiple periods.
# 25.05.2013: Fixed returning original column names for single occurring names.
# 24.05.2013: First version.

#' @title Column Names
#'
#' @description
#' Internal helper function.
#'
#' @details
#' Takes a data frame as input and return either column names
#' occurring once or multiple times. Matching is done by the 'base name'
#' (the substring to the left of the last period, if any). The return type
#' is a string vector by default, or a single string of column names separated
#' by a string 'concatenate' (see 'collapse' in \code{paste} for details).
#' There is an option to limit multiple names to those with a number suffix.
#' 
#' @param data data.frame.
#' @param slim logical, TRUE returns column names occurring once,
#' FALSE returns column names occurring multiple times.
#' @param concatenate string, if not NULL returns a single string with column
#' names concatenated by the provided string instead of a vector.
#' @param numbered logical indicating if repeated column names must have a number suffix.
#' @param debug logical indicating printing debug information.
#' 
#' @return character, vector or string.
#' 
#' @export
#' 
#' @importFrom utils str
#' 

colNames <- function(data, slim=TRUE, concatenate=NULL, numbered=TRUE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("data")
    print(str(data))
    print("slim")
    print(slim)
    print("concatenate")
    print(concatenate)
    print("numbered")
    print(numbered)
  }

  # Initiate variables.
  ret <- NA
  columnNames <- names(data)
  matchingNames <- vector()
  
  if(debug){
    print("columnNames")
    print(columnNames)
  }
  
  if(!is.null(columnNames)){
  
    # Guess base names i.e. "Allele.1" -> "Allele"
    # To the first dot.
    #baseNames <- unique(gsub("(\\w*)\\..*", "\\1", columnNames))
    # To the last dot.
    #baseNames <- unique(gsub("(.*)\\.(.*)", "\\1", columnNames))
    # To the last dot, but return base name for all columns.
    baseNames <- gsub("(.*)\\.(.*)", "\\1", columnNames)
    
    if(debug){
      print("baseNames")
      print(baseNames)
    }
    
    # Create match vector.
    namesCount <- table(baseNames)
    matches <- rep(NA, length(baseNames))
    for(n in seq(along=namesCount)){
      matches[baseNames == names(namesCount[n]) ] <- namesCount[[n]]
    }

    # Indicate column names ending with a number.
    if(numbered){
      number <- suppressWarnings(!is.na(as.numeric(substr(columnNames,
                                                          nchar(columnNames),
                                                          nchar(columnNames)))))
      
      if(debug){
        print("number")
        print(number)
      }
      
    }

    
    if(slim){
      # Return column names matching 'base name' that occured once.
      
      # Get 'base name' that occured once.
      if(numbered){
        singleNames <- unique(baseNames[matches == 1 | !number])
      } else{
        singleNames <- unique(baseNames[matches == 1])
      }

      # Get original column names matching 'base name'.
      for(s in seq(along=singleNames)){
        matchingNames <- c(matchingNames, grep(pattern=singleNames[s], x=columnNames, value = TRUE, fixed = FALSE))
      }
      
      if(is.null(concatenate)){
        # Return as vector.
        ret <- matchingNames
      } else {
        # Return concatenated as single string.
        ret <- paste(matchingNames, collapse=concatenate)
      }
      
    } else {
      # Return column names matching 'base name' that occured multiple times.
      if(is.null(concatenate)){
        # Return as vector.
        if(numbered){
          ret <- unique(baseNames[matches > 1 & number])
        } else {
          ret <- unique(baseNames[matches > 1])
        }
      } else {
        # Return concatenated as single string.
        if(numbered){
          ret <- paste(unique(baseNames[matches > 1 & number]), collapse=concatenate)
        } else {
          ret <- paste(unique(baseNames[matches > 1]), collapse=concatenate)
        }
      }
    }
  }
  
  if(length(ret) == 0){
    ret <- NA
  }
    
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  return(ret)
  
}
