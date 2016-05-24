################################################################################
# TODO LIST
# TODO: option to drop/keep unspecified columns.
# TODO: option to drop/keep 'OL'.

################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added more attributes to result.
# 06.01.2016: Added attributes to result.
# 29.08.2015: Added importFrom.
# 25.08.2015: Fixed error when 'stack' and 'slim' is empty or "".
# 01.06.2015: Fixed columns is found using 'match' instead of 'grep'
#  (fixes problem with partial matching).
# 25.05.2015: Renamed parameters (keepAllFixed -> keep.na)
# 23.05.2014: Improved error message.
# 23.01.2014: Fixed bug when only one column in 'fix'.
# 13.01.2014: Completely re-written for improved performance.
# <13.01.2014: Renamed parameters (slim.col -> stack / fix.col -> fix (as earlier)
# <13.01.2014: Fixed returned factor levels. Added as.matrix to return value.
# <13.01.2014: Renamed parameters (slim -> slim.col / fixed -> fix.col
# <13.01.2014: to avoid function/parameter slim to crash. 
# <13.01.2014: Roxygenized.
# <13.01.2014: new parameter 'keep.na' - WORKING for two key kolumns, but not for unslim
# <13.01.2014: new name flattenGM() -> slim(), new parameter names.

#' @title Slim Data Frames
#'
#' @description
#' Slim data frames with repeated columns.
#'
#' @details
#' Stack repeated columns into single columns.
#' For example, the following data frame:
#'  Sample.Name|Marker|Allele.1|Allele.2|Size.1|Size.2|Data.Point..
#' using this command:
#'  slim(data, fix=c("Sample.Name","Marker"), stack=c("Allele","Size"))
#' would result in this data frame (NB! 'Data.Point' is dropped):
#'  Sample.Name|Marker|Allele|Size
#' 
#' @param data data.frame.
#' @param fix vector of strings with colum names to keep fixed.
#' @param stack vector of strings with colum names to slim.
#' @param keep.na logical, keep a row even if no data.
#' @param debug logical indicating printing debug information.
#' 
#' @return data.frame
#' 
#' @export
#' 
#' @importFrom utils str head
#' 


slim <- function(data, fix=NULL, stack=NULL, 
                 keep.na=TRUE, debug=FALSE){

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("data")
    print(str(data))
    print("fix")
    print(fix)
    print("stack")
    print(stack)
    print("keep.na")
    print(keep.na)
  }
  
  # Check data ----------------------------------------------------------------
  
  # Check parameters.  
  if(!is.logical(keep.na)){
    stop("'keep.na' must be logical",
         call. = TRUE)
  }
  
  if(!is.data.frame(data)){
    stop("'data' must be a data.frame",
         call. = TRUE)
  }

  # Prepare -------------------------------------------------------------------

  # Handle fix="" and fix=character(0)
  if(length(fix) == 0 || all(nchar(fix) == 0)){
    fix = NULL
    message("'fix' set to NULL")
  }
  
  # Handle stack="" and stack=character(0)
  if(length(stack) == 0 || all(nchar(stack) == 0)){
    stack = NULL
    message("'stack' set to NULL")
  }
  
  # Slim ----------------------------------------------------------------------
  
  if(!is.null(stack)){
    
    # Get columns to slim.
    slimCols <- NULL
    for(c in seq(along=stack)){
      slimCols[c] <- list(grepl(stack[c], names(data)))
    }
    
    # Number of columns to slim
    nbSlimCol <- unlist(lapply(slimCols, sum))
    nbCol <- unique(nbSlimCol)
    # Check if equal.
    if(length(nbCol) != 1){
      stop(paste("Columns to stack must have equal number of columns each!",
                 paste(paste(stack, nbSlimCol, sep=":"), collapse="\n"),
                 "The most common problem is multiple columns matching the same 'base' name.",
                 "Here are your column names:",
                 paste(names(data), collapse=", "), sep="\n"))
    }
    
    if(!is.null(fix)){
      
      # Get fixed indexes.
      fixedIndex <- vector()
      for(k in seq(along=fix)){
        #fixedIndex[k] <- grep(fix[k], names(data))
        fixedIndex[k] <- match(fix[k], names(data))
      }
      
      # Get fixed data.
      fixedData <- as.data.frame(data[,fixedIndex])
      names(fixedData) <- fix
      
      # Make a list of matrixes for columns to stack.
      listStack <- list()
      listRep <- vector()
      for(c in seq(along=slimCols)){
        
        # Get columns for current stack key.
        matrixStack <- data[,slimCols[[c]]]
        
        # Count number of values.
        values <- rowSums(!is.na(matrixStack))
        
        # Transpose and vectorize.
        vectorStack <- c(t(matrixStack))
        
        if(keep.na) {
          
          # Keep one value per fixed row.
          values <- replace(values, values==0, 1)
          
          # Create a boolean vector with values to keep.
          bolVec <- rep(c(TRUE,FALSE), nrow(data))
          bolTimes <- vector()
          bolTimes[seq(from=1, to=length(values)*2, by=2)] <- values
          bolTimes[seq(from=2, to=length(values)*2, by=2)] <- nbCol - values
          keep <- rep(bolVec, times=bolTimes)
          
          # Extract values to keep.
          vectorStack <- vectorStack[keep]
          
        } else {
          
          # Extract values to keep.
          vectorStack <- vectorStack[!is.na(vectorStack)]
          
        }
        
        # Add values to list.
        listRep[c] <- list(values)
        
        # Transpose and vectorize matrix, and put in list.
        listStack[c] <- list(vectorStack)
        
      }
      
      # Check if all listRep's are equal.
      for(i in seq(along=listRep)){
        
        # Compare lists.
        if(!identical(listRep[i], listRep[1])){
          
          # Convert mismatch to vectors.
          testA <- unlist(listRep[1])
          testB <- unlist(listRep[i])
          
          # Find row causing the error.
          for(e in seq(along=testA)){
            
            # Compare elements.
            if(testA[e] != testB[e]){
              
              stop(paste("Different repeat patterns detected for stacked columns!\n",
                         "Caused by: ",
                         paste(paste(names(fixedData), ":", fixedData[e,], sep=""),
                               collapse=", ") ,"\n",
                         "Please fix and try again!"), sep="")
              
            }
            
          }
          
        }
      }
      
      # Loop over columns in fixed data.
      fixedDataExt <- list()
      for(k in seq(along=fix)){
        
        # Repeat each 'row' to fit the stacked data.
        fixedDataExt[k] <- list(rep(fixedData[,k], times=unlist(listRep[1])))
        
      }
      
      # Get number of rows.
      numberOfRows <- length(fixedDataExt[[1]])
      # Create a data frame for the result.
      res <- data.frame(matrix(NA, numberOfRows, length(fix) + length(stack)))
      # Add new column names.
      names(res) <- paste(c(fix,stack))
      
      # Combine fixed and stacked list into one.
      listRes <- append(fixedDataExt, listStack)
      for(i in seq(along=listRes)){
        res[,i] <- as.character(listRes[[i]])
      }
      
      # Add attributes to result.
      attr(res, which="slim, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
      attr(res, which="slim, call") <- match.call()
      attr(res, which="slim, date") <- date()
      attr(res, which="slim, data") <- substitute(data)
      attr(res, which="slim, fix") <- fix
      attr(res, which="slim, keep.na") <- keep.na
      attr(res, which="slim, stack") <- stack

      if(debug){
        print(head(res))
        print(str(res))
        print(paste("EXIT:", match.call()[[1]]))
      }
      
      return(res)
      
    } else {

      message("fix=NULL, return data unchanged!")
      
      # Return the data frame unchanged.
      return(data)
      
    }
    
  } else {

    message("stack=NULL, return data unchanged!")
    
    # Return the data frame unchanged.
    return(data)
    
  }
  
}
