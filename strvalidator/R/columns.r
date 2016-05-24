################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 22.12.2015: Added option to add new column by setting col1=NA.
# 28.08.2015: Added importFrom.
# 24.05.2015: First version.

#' @title Column Actions
#'
#' @description
#' Perform actions on columns.
#'
#' @details
#' Perform actions on columns in a data frame. There are five actions:
#' concatenate, add, multiply, subtract, divide. The selected action can be
#' performed on two columns, or one column and a fixed value, or a new column
#' can be added. A target column for the result is specified.
#' NB! if the target column already exist it will be overwritten,
#' else it will be created. A common use is to create a unique 
#' Sample.Name from the existing Sample.Name column and e.g. the File.Name
#' or File.Time columns. It can also be used to calculate the Amount from
#' the Concentration.
#' 
#' @param data a data frame.
#' @param col1 character column name to perform action on.
#' @param col2 character optional second column name to perform action on.
#' @param operator character to indicate operator: & concatenate, + add,
#'  * multiply, - subtract, / divide.
#' @param fixed character or numeric providing the second operand if 'col2'
#'  is not used.
#' @param target character to specify column name for result. Default is to
#' overwrite 'col1'. If not present it will be added.
#' @param debug logical to indicate if debug information should be printed.
#' 
#' @return data frame.
#' 
#' @export
#' 
#' @importFrom utils str
#' 
#' @examples
#' # Get a sample dataset.
#' data(set2)
#' # Add concatenate Sample.Name and Dye.
#' set2 <- columns(data=set2, col1="Sample.Name", col2="Dye")
#' # Multiply Height by 4.
#' set2 <- columns(data=set2, col1="Height", operator="*", fixed=4)
#' # Add a new column.
#' set2 <- columns(data=set2, operator="&", fixed="1234", target="Batch")

columns <- function(data, col1=NA, col2=NA, operator="&",
                    fixed=NA, target=NA, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("data")
    print(str(data))
    print("col1")
    print(col1)
    print("col2")
    print(col2)
    print("operator")
    print(operator)
    print("fixed")
    print(fixed)
    print("target")
    print(target)
  }
  
  # Check data ----------------------------------------------------------------

  if(class(data) != "data.frame"){
    stop("'data' must be of type data.frame",
         call. = TRUE)
  }

  if(!is.na(col1)){
    if(!col1 %in% names(data)){
      stop("'col1' must be a column in 'data'",
           call. = TRUE)
    }
  }  

  if(!is.na(col2)){
    if(!col2 %in% names(data)){
      stop("'col2' must be a column in 'data'",
           call. = TRUE)
    }
  }  
  
  if(!operator %in% c("&", "+", "-", "*", "/")){
    stop("The following operators are supported: &, +, -, *, and /.",
         call. = TRUE)
  }

  # Prepare -------------------------------------------------------------------
  
  # Get values from column.
  if(is.na(col1)){
    # Leave empty to allow adding a new column without specifying a source.
    value1 <- rep("", nrow(data))
    if(operator!="&"){
      operator <- "&"
      message("'operator' changed to '&', since column 1 is not specified.")
    }
  } else {
    # Get source values from column 1.
    value1 <- data[,col1]
  }

  # Get second values from column or argument.
  if(col2 %in% names(data)){
    value2 <- data[,col2]
  } else {
    value2 <- fixed
  }

  # If no target overwrite column 1.
  if(is.na(target)){
    target <- col1
  }

  # Perform actions -----------------------------------------------------------
  
  # Select action.
  if(operator == "&"){
    
    valueNew <- paste(value1, value2, sep="")
    
  } else if(operator == "+"){
    
    value1 <- as.numeric(value1)
    value2 <- as.numeric(value2)
    valueNew <- value1 + value2
    
  } else if(operator == "*"){
    
    value1 <- as.numeric(value1)
    value2 <- as.numeric(value2)
    valueNew <- value1 * value2
    
  } else if(operator == "-"){
    
    value1 <- as.numeric(value1)
    value2 <- as.numeric(value2)
    valueNew <- value1 - value2
    
  } else if(operator == "/"){
    
    value1 <- as.numeric(value1)
    value2 <- as.numeric(value2)
    valueNew <- value1 / value2
    
  }
  
  # Add new values to data frame.
  data[target] <- valueNew
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result.
  return(data)
  
}
