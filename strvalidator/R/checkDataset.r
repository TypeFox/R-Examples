################################################################################
# TODO LIST
# TODO: Change parameter names to format: lower.case

################################################################################
# CHANGE LOG (last 20 changes)
# 04.05.2015: 'slimcol' and 'stringcol' now accept vectors. 
# 06.05.2014: First version.

#' @title Check Dataset
#'
#' @description
#' Internal function to check a data.frame before analysis.
#'
#' @details Check that the object exist, there are rows, the required columns exist,
#' if data.frame is 'fat', and if invalid strings exist. Show error message if not.
#' 
#' @param name character name of data.frame.
#' @param reqcol character vector with required column names.
#' @param slim logical TRUE to check if 'slim' data.
#' @param slimcol character vector with column names to check if 'slim' data.
#' @param string character vector with invalid strings in 'stringcol', return FALSE if found.
#' @param stringcol character vector with column names to check for 'string'.
#' @param env environment where to look for the data frame.
#' @param parent parent gWidget.
#' @param debug logical indicating printing debug information.
#' 

checkDataset <- function(name, reqcol=NULL, slim=FALSE, slimcol=NULL,
                         string=NULL, stringcol=NULL,
                         env=parent.frame(), parent=NULL, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("name")
    print(name)
    print("reqcol")
    print(reqcol)
    print("slim")
    print(slim)
    print("stringcol")
    print(stringcol)
    print("env")
    print(environmentName(env))
    print("parent")
    print(parent)
  }
  
  # Initiate variables.
  ok <- TRUE
  messageText <- NULL
  
  # Check if dataset exist in environment.
  if(exists(name, envir=env, inherits = FALSE)){
    
    # Get dataset.
    df <- get(name, envir=env)

    # Check if empty dataset.
    if(nrow(df) == 0){
      
      # Construct error message.
      messageText <- c("Dataset contain no rows!")
      
      # Change flag.
      ok <- FALSE
      
    } else if(!is.null(reqcol) & !all(reqcol %in% colnames(df))){
      # Check for required column names.
      
      missingCol <- reqcol[!reqcol %in% colnames(df)]
      
      # Construct error message.
      messageText <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      # Change flag.
      ok <- FALSE
    
    } else if(slim & !is.null(slimcol)){
      
      # Loop over columns to check.
      for(c in seq(along=slimcol)){

        # Check if slimmed.
        slimmed <- sum(grepl(slimcol[c], names(df), fixed=TRUE)) == 1
        
        if(!slimmed){
          
          # Construct error message.
          messageText <- paste("The dataset is too fat!\n\n",
                               "There can only be 1", slimcol[c], "column\n",
                               "Slim the dataset", sep="")
          
          # Change flag.
          ok <- FALSE
          
        }
        
      }
      
    } else if(!is.null(string) & !is.null(stringcol)){
      
      # Loop over columns to check.
      for(c in seq(along=stringcol)){

        if(any(string %in% df[,stringcol[c]])){
          
          # Construct error message.
          messageText <- paste("'", string, "' detected in column ", stringcol[c], "!\n",
                               "Please make sure that data is clean/filtered.\n",
                               sep="")
          
          # Change flag.
          ok <- FALSE
          
        }
        
      }
      
    } else {
      
      # Dataset passed!
      ok <- TRUE
      
    }
    
  } else {
 
    # Construct error message.
    messageText <- NULL # NB! Can't show error message because "<Select dataset>" is in drop-downs.  
    
    # Change flag.
    ok <- FALSE
    
  }

  if(debug){
    print("ok")
    print(ok)
    print("messageText")
    print(messageText)
  }
  
  # Show error message.
  if(!ok & !is.null(messageText)){
    if(is.null(parent)){
      # Write in command prompt.
      message(messageText)
    } else {
      # Show message box.
      gmessage(messageText, title="message", icon = "error", parent = parent) 
    }
  }
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result.
  return(ok)
  
}